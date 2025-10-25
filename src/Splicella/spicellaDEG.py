import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import itertools


class SplicellaDEG_sampling:
    def __init__(self, genecount_anndata, 
                by='celltype', layer=None, pairwise=True,
                gene_name='gene_name', 
                 method='wilcoxon',
                 valid_gene_thr=5, sample_size = 900, repeat = 1000, 
                 pval_column = 'pvals_adj', pval_cutoff=0.05, lfc_cutoff=1,
                 seed=42):
        self.gad = genecount_anndata.copy()
        self.gene_list = self.gad.var_names
        
        self.state = by
        self.states = self.gad.obs[self.state].cat.categories
        self.pairwise = pairwise
        self.layer = layer
        self.gene_name = gene_name
        self.method = method
        
        self.valid_gene_thr = valid_gene_thr
        self.sample_size = sample_size
        self.repeat = repeat

        self.pval_column = pval_column
        self.pval_cutoff = pval_cutoff
        self.lfc_cutoff = lfc_cutoff

        np.random.seed(seed)

        if pairwise:
            self.gene_detect_df = self.get_gene_detect_df(self.state)
        else:
            self.prepare_nonpairwise()

    def prepare_nonpairwise(self):
        for state in self.states:
            self.gad.obs[f'celltype_{state}'] = [state if bools else 'others' for bools in self.gad.obs[self.state]==state]
            
        
    def get_gene_detect_df(self, by):
        nonzero_gene_adata = sc.get.aggregate(self.gad, by=by,
                                              layer=self.layer, func='count_nonzero')
    
        gene_detect_df = pd.DataFrame(nonzero_gene_adata.layers['count_nonzero'], 
                                      index = nonzero_gene_adata.obs_names, 
                                      columns=nonzero_gene_adata.var_names)
        
        # percentage of cells expressing the gene
        cell_counts = self.gad.obs[self.state].value_counts().loc[gene_detect_df.index]
        gene_detect_df = gene_detect_df.T * 100 / cell_counts.values
        return gene_detect_df
        
    
    def get_deg_df(self, adata, by, state1, state2, 
                   gene_detect_df, pval_cutoff=None, lfc_cutoff=None):
        if pval_cutoff is None:
            pval_cutoff = self.pval_cutoff
        if lfc_cutoff is None:
            lfc_cutoff = self.lfc_cutoff
        results = []
        state1_genes = gene_detect_df.loc[:, state1] >= self.valid_gene_thr
        state2_genes = gene_detect_df.loc[:, state2] >= self.valid_gene_thr
        
        filtered_genes = gene_detect_df[state1_genes | state2_genes].index

        filtered_adata = self.gad[:, self.gad.var[self.gene_name].isin(filtered_genes)].copy()

        
        sc.tl.rank_genes_groups(filtered_adata, 
                                by,
                                groups=[state1], reference=state2, 
                                method=self.method, layer=self.layer)
        res = sc.get.rank_genes_groups_df(filtered_adata, group=state1)
        res['compare'] = f'{state1}_{state2}'
        res['DEG_UP'] = (res[self.pval_column]<=pval_cutoff)&(res['logfoldchanges']>=lfc_cutoff)
        return res


    def sampled_deg_computation(self, by, state1, state2):
        deg_df_all = []
        
        if self.pairwise:
            gene_detect_df = self.gene_detect_df
        else:
            gene_detect_df = self.get_gene_detect_df(by=f'celltype_{state1}')


        for i in range(self.repeat):
            total_size = [self.gad[self.gad.obs[by]==state1].shape[0], 
                          self.gad[self.gad.obs[by]==state2].shape[0]]
            minimum_size = min(total_size)
            sample_size = min(int(minimum_size*0.9), self.sample_size)
   
            state1_ids = self.gad.obs[self.gad.obs[by]==state1].sample(n=sample_size).index
            state2_ids = self.gad.obs[self.gad.obs[by]==state2].sample(n=sample_size).index
            target_ids = state1_ids.union(state2_ids)
                
            subset_adata = self.gad[target_ids, :].copy()

            deg = self.get_deg_df(subset_adata, by, state1, state2, 
                                  gene_detect_df)
            deg["r"] = i
            deg_df_all.append(deg)

        deg_df_all = pd.concat(deg_df_all, ignore_index=True)
        return deg_df_all, gene_detect_df

    def DEA_result_per_gene(self, deg_df_all):

        DEG_df = deg_df_all.pivot_table(columns='r', index='names', values='DEG_UP')
        lfc_df = deg_df_all.pivot_table(columns='r', index='names', values='logfoldchanges')
        pval_df = deg_df_all.pivot_table(columns='r', index='names', values='pvals_adj')


        return DEG_df, lfc_df, pval_df

    def compute_tau(self, gene_detect_df, DEG_df, by, state1, state2):
        valid_genes = gene_detect_df[(gene_detect_df > self.valid_gene_thr).any(axis=1)].index

        nosampled_DEG = self.get_deg_df(self.gad, by, state1, state2, gene_detect_df)

        PFER = len(nosampled_DEG[nosampled_DEG['DEG_UP']])
        
        Lrs = DEG_df.sum(axis=0)
        q = np.mean(Lrs)
        p = len(valid_genes)
        tau = 0.5 * (q**2 / (PFER * p) + 1)
        
        return(Lrs, tau, q, PFER, nosampled_DEG)

    def robustDEGs(self, gene_detect_df, nosampled_DEG, DEG_df, pval_df, lfc_df):
        valid_genes = gene_detect_df[(gene_detect_df > self.valid_gene_thr).any(axis=1)].index
        nosampled_DEG_ = nosampled_DEG[nosampled_DEG['DEG_UP']]

        Pi_mean = DEG_df.mean(axis=1)
        pval_mean = pval_df.mean(axis=1)
        lfc_mean = lfc_df.mean(axis=1)

        ups = Pi_mean[Pi_mean>tau].index

        up_diff = set(nosampled_DEG['names']).difference(ups)
        up_common = set(nosampled_DEG['names']).intersection(ups)

        up_diff_dic = {}
        up_common_dic = {}

        records = []

        for _, line in nosampled_DEG_.iterrows():
            name, pval, logfc = line['names'], line[self.pval_column], line['logfoldchanges']
            if (ups==name).any():
                set_ = 'common'
            else:
                set_ = 'diff'
            records.append([name, pval, logfc, set_])
        robustDEG_summary = pd.DataFrame(records, columns=['name', self.pval_column, 'logfoldchnages', 'set'])

        return robustDEG_summary

    def splicelladeg(self, groups=None, ):
        if groups is not None:
            self.states = list(set(groups)&set(self.states))

        deg_df_alls = []
        Lrs_dict, tau_dict, q_dict, pfer_dict = {}, {}, {}, {}
        robustDEGs = []
        if self.pairwise:
            for state, state_conv in itertools.permutations(self.states, 2):
                deg_df_all, gene_detect_df = self.sample_deg_computation(self.state, state, state_conv)
                DEG_df, lfc_df, pval_df = self.DEA_result_per_gene(deg_df_all)
                Lrs, tau, q, PFER, nosampled_DEG = self.compute_tau(gene_detect_df, DEG_df, self.state, state, state_conv)
                robustDEG_summary = self.robustDEGs(gene_detect_df, nosampled_DEG, DEG_df, pval_df, lfc_df)    

                deg_df_alls.append(deg_df_all)
                compare = f'{state}_{state_conv}'
                Lrs_dict[compare] = Lrs
                tau_dict[compare] = tau
                q_dict[compare] = q
                PFER_dict[compare] = PFER
                robustDEGs.append(robustDEG_summary)
                
                
        else:
            for state in self.states:
                state_conv = 'others'
                by = f'celltype_{state}'
                deg_df_all, gene_detect_df = self.sample_deg_computation(self.state, state, state_conv)
                DEG_df, lfc_df, pval_df = self.DEA_result_per_gene(deg_df_all)
                Lrs, tau, q, PFER, nosampled_DEG = self.compute_tau(gene_detect_df, DEG_df, by, state, state_conv)
                robustDEG_summary = self.robustDEGs(gene_detect_df, nosampled_DEG, DEG_df, pval_df, lfc_df)    

                deg_df_alls.append(deg_df_all)
                compare = f'{state}_{state_conv}'
                Lrs_dict[compare] = Lrs
                tau_dict[compare] = tau
                q_dict[compare] = q
                PFER_dict[compare] = PFER
                robustDEGs.append(robustDEG_summary)
        deg_df_alls = pd.concat(deg_df_alls)
        robustDEGs = pd.concat(robustDEGs)
        return(deg_df_alls, robustDEGs, Lrs_dict, tau_dict, q_dict, PFER_dict)
