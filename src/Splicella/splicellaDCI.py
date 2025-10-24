import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
import itertools

class SplicellaDCI_pairwise:
    def __init__(self,
                isoform_anndata, 
                by='celltype', 
                layer='total_raw',
                gene_name='gene_name',
                ):
        self.iad = isoform_anndata
        self.state = by
        self.states = self.iad.obs[self.state].cat.categories
        self.layer = layer
        self.gene_name = gene_name
        self.gene_iso_num = pd.DataFrame(np.unique(self.iad.var[self.gene_name], return_counts=True)).T
        self.multiple_gene = self.gene_iso_num[self.gene_iso_num[1]>=2][0]
        self.iad = self.iad[:, self.iad.var[self.gene_name].isin(self.multiple_gene)]


    def get_isoform_preval(self):
        """
        get calcuated prevalance and raw count per cell/state
        """
        
        self.iso_agg = sc.get.aggregate(self.iad, by=self.state, layer=self.layer, func='sum')
        
        gene_agg = sc.get.aggregate(self.iad, by=self.gene_name, axis=1, layer=self.layer, func='sum')
        gene_agg = sc.get.aggregate(gene_agg, by=self.state, layer='sum', func='sum')
    
        genecount_agg = sc.get.aggregate(self.iad, by=self.gene_name, layer=self.layer,  func='count_nonzero', axis=1)
        genecount_agg = sc.get.aggregate(genecount_agg, by=self.state, layer='count_nonzero', func='count_nonzero')
    
        self.gene_sum_df = pd.DataFrame(gene_agg.layers['sum'], 
                               index=gene_agg.obs_names,
                               columns=gene_agg.var_names)
        self.gene_detect_df = pd.DataFrame(genecount_agg.layers['count_nonzero'],
                                  index=genecount_agg.obs_names,
                                  columns=genecount_agg.var_names)
    
        self.gene_detect_percent_df = self.gene_detect_df.T*100/np.unique(self.iad.obs[self.state],
                                                                  return_counts=True)[1]
        
        total_gc_forprev = []
        for _, line in self.iad.var.iterrows():
            gene_name = line[self.gene_name]
            gene_sum = self.gene_sum_df[gene_name]
            gene_sum[gene_sum==0] = 1
            total_gc_forprev.append(gene_sum.values.ravel())
        self.iso_agg.layers['preval'] = self.iso_agg.layers['sum']/np.array(total_gc_forprev).T
        
        return None

    def DCI_test(self):
        valid_genes = self.gene_detect_percent_df[
        (self.gene_detect_percent_df >= self.vg_threshold).any(axis=1)].index
        chitest_results = []
        isoform_preval_result = []
        
        for state, state_conv in itertools.permutations(self.states, 2):
            
            state_genes = self.gene_detect_percent_df.loc[valid_genes, state] >= self.vg_threshold
            state_conv_genes = self.gene_detect_percent_df.loc[valid_genes, state_conv] >= self.vg_threshold
    
            filtered_genes = valid_genes[state_genes & state_conv_genes]
            chitest_result = []
            for gene in filtered_genes:
                genesum = self.gene_sum_df[gene]
                isoforms = self.iso_agg.var[self.iso_agg.var[self.gene_name]==gene].index
                state_isoforms = self.iso_agg[state, isoforms].layers['sum']
                state_conv_isoforms = self.iso_agg[state_conv, isoforms].layers['sum']
    
                observed = np.concatenate([state_isoforms, state_conv_isoforms], axis=0)
                observed_sum = observed.sum(axis=0)
                
                observed = observed[:, observed_sum>=self.vi_threshold]
                isoforms = isoforms[observed_sum>=self.vi_threshold]
    
                if len(isoforms)<2:
                    continue
                    
                try:
                    chi2_stat, p_val, _, _ = chi2_contingency(observed)
                except:
                    continue

                state_preval = self.iso_agg[state, isoforms].layers['preval'].ravel()*100
                state_conv_preval = self.iso_agg[state_conv, isoforms].layers['preval'].ravel()*100
                preval_diff = state_preval - state_conv_preval

                
                chitest_result.append({
                    'gene': gene,
                    'state': state,
                    'state_conv': state_conv,
                    'chi2_statistic': chi2_stat,
                    'p_value': p_val,
                    'max_prevaldiff': preval_diff.max(),
                    'max_abs_prevaldiff': abs(preval_diff).max()
                })
    
                for i, (isoform, p_d) in enumerate(zip(isoforms, preval_diff)):
                    isoform_preval_result.append({
                        'isoform': isoform,
                        'gene': gene,
                        'state': state,
                        'state_conv': state_conv,
                        'state_preval': state_preval[i],
                        'conv_preval': state_conv_preval[i],
                        'prevaldiff': p_d,
                    }) 
            chitest_result = pd.DataFrame(chitest_result)
            rejected, pvals_corrected, _, _ = multipletests(chitest_result['p_value'],
                                                           method='fdr_bh')
            chitest_result['FDR'] = pvals_corrected
            chitest_results.append(chitest_result)
                    
        self.chitest_df = pd.concat(chitest_results)
        self.isoform_preval_df = pd.DataFrame(isoform_preval_result)

        
    def get_dcidf(self,
                  pval_cutoff=0.05, pval_column='FDR', prevaldiff_cutoff=10, use_absolute_threshold=False):
        self.pval_cutoff = pval_cutoff
        self.pval_column = pval_column
        self.prevaldiff_cutoff = prevaldiff_cutoff
        self.use_absolute_threshold = use_absolute_threshold

        if self.use_absolute_threshold:
            self.DCI_df = self.chitest_df[
            (self.chitest_df[self.pval_column]<=self.pval_cutoff) &\
            (self.chitest_df['max_abs_prevaldiff']>=self.prevaldiff_cutoff)]


        else:
            self.DCI_df = self.chitest_df[
            (self.chitest_df[self.pval_column]<=self.pval_cutoff) &\
            (self.chitest_df['max_prevaldiff']>=self.prevaldiff_cutoff)]

        DCI_iso_df = []
        for _, line in self.DCI_df.groupby(['gene', 'state']):
            gene, state = _
            isodf = self.isoform_preval_df[(self.isoform_preval_df['gene']==gene)&(self.isoform_preval_df['state']==state)]
            isodf = isodf[isodf['prevaldiff']>=10]
            DCI_iso_df.append(isodf)
        self.DCI_iso_df = pd.concat(DCI_iso_df)

    def splicelladci(self, 
                     valid_gene_threshold=5, 
                     valid_isoform_threshold=10):
        self.vg_threshold = valid_gene_threshold
        self.vi_threshold = valid_isoform_threshold
        
        self.get_isoform_preval()
        self.DCI_test()


class SplicellaDCI_nonpairwise:
    def __init__(self,
                isoform_anndata, 
                by='celltype', 
                layer='total_raw',
                gene_name='gene_name'):
        self.iad = isoform_anndata
        self.state = by
        self.states = self.iad.obs[by].cat.categories
        self.layer = layer
        self.gene_name = gene_name

        self.gene_iso_num = pd.DataFrame(np.unique(self.iad.var[self.gene_name], return_counts=True)).T
        self.multiple_gene = self.gene_iso_num[self.gene_iso_num[1]>=2][0]
        self.iad = self.iad[:, self.iad.var[self.gene_name].isin(self.multiple_gene)]

    def get_isoform_preval(self):
        """
        get calcuated prevalance and raw count per cell/state
        """
        #isoform count
        self.iso_agg = sc.get.aggregate(self.iad, by=self.state, layer=self.layer, func='sum') #isoformcount by celltype
        comp_isoagg = np.zeros(self.iso_agg.shape) #isoformcout by celltype for compare
        for i in range(self.iso_agg.n_obs):
            agg_comp_i = np.delete(self.iso_agg.layers['sum'], i, axis=0)
            comp_isoagg[i] = np.sum(agg_comp_i, axis=0)
        self.iso_agg.layers['sum_comp'] = comp_isoagg

        #genecount
        self.gene_agg = sc.get.aggregate(self.iad, by=self.gene_name, axis=1, layer=self.layer, func='sum')
        self.gene_agg = sc.get.aggregate(self.gene_agg, by=self.state, layer='sum', func='sum') #genecount by celltype

        comp_geneagg = np.zeros(self.gene_agg.shape)
        for i in range(self.gene_agg.n_obs):
            agg_comp_i = np.delete(self.gene_agg.layers['sum'], i, axis=0)
            comp_geneagg[i] = np.sum(agg_comp_i, axis=0)
        self.gene_agg.layers['sum_comp'] = comp_geneagg


        #percent of gene detection per cell
        genecount_agg = sc.get.aggregate(self.iad, by=self.gene_name, layer=self.layer,  func='count_nonzero', axis=1)
        genecount_agg = sc.get.aggregate(genecount_agg, by=self.state, layer='count_nonzero', func='count_nonzero') #gene exist by celltype

        self.gene_detect_df = pd.DataFrame(genecount_agg.layers['count_nonzero'],
                                  index=genecount_agg.obs_names,
                                  columns=genecount_agg.var_names)
        comp_count_nonzero = np.zeros_like(genecount_agg.layers['count_nonzero'])
        for i in range(genecount_agg.n_obs):
            agg_without_i = np.delete(genecount_agg.layers['count_nonzero'], i, axis=0)
            comp_count_nonzero[i] = np.sum(agg_without_i, axis=0)
            
        self.gene_detect_comp_df = pd.DataFrame(comp_count_nonzero,
                                  index=genecount_agg.obs_names,
                                  columns=genecount_agg.var_names)
        
        self.celltype_counts = self.iad.obs[self.state].value_counts()[self.gene_detect_df.index]
        self.celltype_comp_counts = self.iad.n_obs - self.iad.obs[self.state].value_counts()[self.gene_detect_comp_df.index]

        self.gene_detect_percent_df = self.gene_detect_df.T*100/self.celltype_counts[genecount_agg.obs_names]
        self.gene_detect_comp_percent_df = self.gene_detect_comp_df.T*100/self.celltype_comp_counts[genecount_agg.obs_names]

        #isoform prevalance
        total_gc_forprev = []
        total_compgc_forprev = []
        for _, line in self.iad.var.iterrows():
            gene_name = line[self.gene_name]
            
            gene_sum = self.gene_agg[:, gene_name].layers['sum'].ravel().copy()
            gene_sum[gene_sum==0] = 1
            total_gc_forprev.append(gene_sum)

            gene_compsum = self.gene_agg[:, gene_name].layers['sum_comp'].ravel().copy()
            gene_compsum[gene_compsum==0] = 1
            total_compgc_forprev.append(gene_compsum)
        self.iso_agg.layers['preval'] = self.iso_agg.layers['sum']/np.array(total_gc_forprev).T
        self.iso_agg.layers['preval_comp'] = self.iso_agg.layers['sum_comp']/np.array(total_compgc_forprev).T

    def DCI_test(self):
        valid_genes = self.gene_detect_percent_df[
        (self.gene_detect_percent_df >= self.vg_threshold).any(axis=1)].index

        chitest_results = []
        isoform_preval_result = []
        
        for state in self.states:
            state_genes = self.gene_detect_percent_df.loc[valid_genes, state] >= self.vg_threshold
            nonstate_genes = self.gene_detect_comp_percent_df.loc[valid_genes, state] >= self.vg_threshold
            filtered_genes = valid_genes[state_genes&nonstate_genes]

            chitest_result = []
            for gene in filtered_genes:
                isoforms = self.iso_agg.var[self.iso_agg.var[self.gene_name]==gene].index
                state_isoforms = self.iso_agg[state, isoforms].layers['sum']
                compstate_isoforms = self.iso_agg[state, isoforms].layers['sum_comp']

                observed = np.concatenate([state_isoforms, compstate_isoforms], axis=0)
                observed_sum = observed.sum(axis=0)
                observed = observed[:, observed_sum>=self.vi_threshold]
                isoforms = isoforms[observed_sum>=self.vi_threshold]

                if len(isoforms)<2:
                    continue

                try:
                    chi2_stat, p_val, _, _ = chi2_contingency(observed)
                except:
                    continue

                state_preval = self.iso_agg[state, isoforms].layers['preval'].ravel()*100
                state_comp_preval = self.iso_agg[state, isoforms].layers['preval_comp'].ravel()*100
                preval_diff = state_preval - state_comp_preval

                chitest_result.append({
                    'gene': gene,
                    'state': state,
                    'chi2_statistic': chi2_stat,
                    'p_value': p_val,
                    'max_prevaldiff': preval_diff.max(),
                    'max_abs_prevaldiff': abs(preval_diff).max()
                })

                for i, (isoform, p_d) in enumerate(zip(isoforms, preval_diff)):
                    isoform_preval_result.append({
                        'isoform': isoform,
                        'gene': gene,
                        'state': state,
                        'state_preval': state_preval[i],
                        'conv_preval': state_comp_preval[i],
                        'prevaldiff': p_d,
                        'abs_prevaldiff': abs(p_d),
                    }) 
            chitest_result = pd.DataFrame(chitest_result)
            rejected, pvals_corrected, _, _ = multipletests(chitest_result['p_value'],
                                                           method='fdr_bh')
            chitest_result['FDR'] = pvals_corrected
            chitest_results.append(chitest_result)
                    
        self.chitest_df = pd.concat(chitest_results)
        self.isoform_preval_df = pd.DataFrame(isoform_preval_result)

        
    def get_dcidf(self,
                 pval_cutoff=0.05, pval_column='FDR', prevaldiff_cutoff=10, use_absolute_threshold=False):
        self.pval_cutoff = pval_cutoff
        self.pval_column = pval_column
        self.prevaldiff_cutoff = prevaldiff_cutoff
        self.use_absolute_threshold = use_absolute_threshold
        if self.use_absolute_threshold:
            self.DCI_df = self.chitest_df[
            (self.chitest_df[self.pval_column]<=self.pval_cutoff) &\
            (self.chitest_df['max_abs_prevaldiff']>=self.prevaldiff_cutoff)]

        else:
            self.DCI_df = self.chitest_df[
            (self.chitest_df[self.pval_column]<=self.pval_cutoff) &\
            (self.chitest_df['max_prevaldiff']>=self.prevaldiff_cutoff)]

        DCI_iso_df = []
        for _, line in self.DCI_df.groupby(['gene', 'state']):
            gene, state = _
            isodf = self.isoform_preval_df[(self.isoform_preval_df['gene']==gene)&(self.isoform_preval_df['state']==state)]
            isodf = isodf[isodf['prevaldiff']>=10]
            DCI_iso_df.append(isodf)
        self.DCI_iso_df = pd.concat(DCI_iso_df)        

    def splicelladci(self, 
                     valid_gene_threshold=5, 
                     valid_isoform_threshold=10):
                     
        self.vg_threshold = valid_gene_threshold
        self.vi_threshold = valid_isoform_threshold        
        self.get_isoform_preval()
        self.DCI_test()



        
        

    




