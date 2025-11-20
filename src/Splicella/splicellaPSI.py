import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.metrics.pairwise import nan_euclidean_distances
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, optimal_leaf_ordering, dendrogram
from scipy.cluster.hierarchy import leaves_list
import itertools
import seaborn as sns

class SplicellaPSI:
    def __init__(self,
                exon_anndata, by,
                include_layer='include_raw', total_layer='total_raw', 
                pseudocount=0.5, mintotalsum=10):
        self.ead = exon_anndata
        self.state = by
        self.include_layer = include_layer
        self.total_layer = total_layer
        self.pseudocount = pseudocount
        self.min_total_per_group = mintotalsum

        if not hasattr(self.ead, 'obs') or self.state not in self.ead.obs.columns:
            raise ValueError(f"The provided by column {self.state} does not exist in the adata object.")

    def compute_group_psi(self):
        """
        get PSI per exon
        """
        agg_inc = sc.get.aggregate(self.ead, by=self.state, layer=self.include_layer, func='sum')
        agg_tot = sc.get.aggregate(self.ead, by=self.state, layer=self.total_layer,   func='sum')
    
        inc = agg_inc.layers['sum']
        tot = agg_tot.layers['sum']
    
        with np.errstate(divide='ignore', invalid='ignore'):
            psi = (inc + self.pseudocount) / (tot + 2*self.pseudocount)
        psi[(tot < self.min_total_per_group)] = np.nan
    
        self.psi_df  = pd.DataFrame(psi, index=agg_inc.obs_names, columns=agg_inc.var_names)
        self.tot_df  = pd.DataFrame(tot, index=agg_inc.obs_names, columns=agg_inc.var_names)

    def compute_all_pairwise_delta_psi(self):
        """
        compute d(PSI) of all group from psi_df (groups × exons)
        """
        self.deltas_dict = {}
        groups = list(self.psi_df.index)
        for g1, g2 in itertools.combinations(groups, 2):
            d = (self.psi_df.loc[g1] - self.psi_df.loc[g2]).abs()
            self.deltas_dict[(g1, g2)] = d

    def call_hVE_from_deltas(self, delta_thr=.5):
        """
        set as hVE if any comparison shows higher ΔPSI than thershold
        return:
          hve_flag: pd.Series(bool, index=exons)   #wheter hVE
          max_delta: pd.Series(float, index=exons)  # maximum ΔPSI for each exon
          n_compared: pd.Series(int, index=exons)   # the number of compair except NaN
        """
        self.compute_group_psi()
        self.compute_all_pairwise_delta_psi()
        
        all_exons = None
        for _k, s in self.deltas_dict.items():
            all_exons = s.index if all_exons is None else all_exons.union(s.index)
        all_exons = pd.Index(sorted(all_exons))
    
        max_delta = pd.Series(0.0, index=all_exons)
        n_compared = pd.Series(0, index=all_exons)
    
        for (_g1,_g2), s in self.deltas_dict.items():
            valid = s.dropna()
            n_compared[valid.index] += 1
            to_update = valid.index[valid.values > max_delta.loc[valid.index].values]
            max_delta.loc[to_update] = valid.loc[to_update].values
    
        hve_flag = (max_delta >= delta_thr) & (n_compared >= 1) 
        return hve_flag, max_delta, n_compared


def get_pseudogroup(adata, time_key='velocity_pseudotime', n_groups=5,
                   replace=True):
    pseudotime = adata.obs[time_key]
    p_group = pd.qcut(pseudotime, q=n_groups, labels=False)
    return(p_group)

def get_optimalgroup(adata, time_key='velocity_pseudotime', group_min=3, group_max=10,
                     include_layer='include_raw', total_layer='total_raw', 
                pseudocount=0.5, mintotalsum=10, delta_thr=.5):
    for i in range(group_min, group_max+1):
        adata.obs[f'p_group{i}'] = get_pseudogroup(adata, time_key, i)

    results_by_k = {}
    metrics_by_k = []
    hve_sets_by_k = {}
    
    for k in range(group_min, group_max+1):
        group_col = f"p_group{k}"
        spcPSI = SplicellaPSI(adata, group_col,
                           include_layer=include_layer, total_layer=total_layer,
                           pseudocount=pseudocount, mintotalsum=mintotalsum)
        hve_flag, max_delta, n_compared = spcPSI.call_hVE_from_deltas(delta_thr)
    
        res = pd.DataFrame({
            'max_delta_psi': max_delta,
            'n_compared_pairs': n_compared,
            'hVE': hve_flag
        })
        results_by_k[k] = {
            'psi': spcPSI.psi_df,
            'tot': spcPSI.tot_df,
            'deltas': spcPSI.deltas_dict,
            'summary': res
        }
    
        hve_sets_by_k[k] = set(res.index[hve_flag])
    
        hve = pd.DataFrame(spcPSI.deltas_dict)[hve_flag]
        min_obs_ratio_cols = 0.6
        
        col_keep = hve.notna().mean(axis=0) >= min_obs_ratio_cols
    
        metrics_ = {'keepcol': hve.loc[:, col_keep].shape[1]/hve.shape[1]}
        metrics_by_k.append(metrics_)
        
    k_metrics_df = pd.DataFrame(metrics_by_k, index=range(group_min, group_max+1))

    return results_by_k, hve_sets_by_k, k_metrics_df

def get_trimmed_eve(df, min_obs_ratio_rows=.6, min_obs_ratio_cols=.6):
    row_keep = df.notna().mean(axis=1) >= min_obs_ratio_rows
    col_keep = df.notna().mean(axis=0) >= min_obs_ratio_cols
    return df.loc[row_keep, col_keep]


def get_linkage(psi_df, hVE_deltas, linkagemethod='ward'):
    psi_df = psi_df.T
    psi_df_ = psi_df.loc[hVE_deltas.index]
    Distance = nan_euclidean_distances(psi_df_.values)
    psi_df_ = psi_df_.iloc[np.setdiff1d(np.arange(psi_df_.shape[0]),
                          np.unique(np.argwhere(np.isnan(Distance)).ravel()))]
    Distance = nan_euclidean_distances(psi_df_.values)
    Z = linkage(squareform(Distance, checks=False), method=linkagemethod)
    return Z
    

def get_cluster(psi_df, hVE_deltas,
                Z, k, clustercriterion='maxclust', clustermap=False,
               ):
    psi_df_ = psi_df.T
    psi_df_ = psi_df_.loc[hVE_deltas.index]
    labels_exon = fcluster(Z, t=k, criterion=clustercriterion)
    exon_labels = pd.Series(labels_exon, index=psi_df_.index, name='cluster')
    if clustermap:
        leaf_order = leaves_list(Z)
        ordered = psi_df_.iloc[leaf_order]
        unique_clusters = sorted(exon_labels.unique())
        palette = sns.color_palette("Set2", n_colors=len(unique_clusters))
        lut = dict(zip(unique_clusters, palette))
        
        row_colors = exon_labels.loc[ordered.index].map(lut)
        
        cg = sns.clustermap(
            ordered.fillna(0), 
            cmap="Spectral_r", 
            row_colors=row_colors, 
            row_cluster=False, 
            col_cluster=False,
            mask=ordered.isna(),
            linewidths=0,
            figsize=(12, 10),
            cbar_kws=dict(label='PSI'),
            yticklabels=False,
        )
        
        return exon_labels, cg
    else:
        return exon_labels
    

