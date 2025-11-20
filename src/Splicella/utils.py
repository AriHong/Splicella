import pandas as pd
import scanpy as sc
import numpy as np
import pyranges as pr
import anndata

def get_annotr_fromgtf(gtf_path, save_path=None):
    anno = pr.read_gtf('refs/gencode.v48.primary_assembly.annotation.gtf.gz').df
    anno_tr = anno[anno['Feature'] == 'transcript'].copy()
    
    anno_exon = anno[anno['Feature'] == 'exon'].copy()
    anno_exon.loc[:, 'Length'] = anno_exon['End'] - anno_exon['Start'] + 1
    tr_length = anno_exon.pivot_table(index='transcript_id', values='Length', aggfunc='sum')
    tr_exonnum = anno_exon.pivot_table(index='transcript_id', values='Length', aggfunc='count')
    
    anno_tr_ = anno_tr[['transcript_name', 'transcript_id', 'gene_name', 'gene_id',
                    'gene_type', 'transcript_type']]

    anno_tr_ = anno_tr_.set_index('transcript_id')
    anno_tr_['length'] = tr_length['Length']
    anno_tr_['exon_num'] = tr_exonnum['Length']
    
    if save_path is not None:
        anno_tr_.to_csv(save_path, index=None)
    return anno_tr_
    
def get_adata_from_iso(file_path, gtf=None, anno_tr=None, save_path=None):
    if anno_tr is None:
        if gtf is None:
            raise ValueError("")
        else:
            anno_tr = get_annotr_fromgtf(gtf)
    else:
        if anno_tr.index.name == 'transcript_id':
            anno_tr = anno_tr.reset_index()

    df = pd.read_csv(file_path, sep='\t')

    var = df['feature_id']
    var = pd.DataFrame(var.values, columns=['transcript_id']).set_index('transcript_id')

    var = pd.merge(var, anno_tr, left_index=True, right_on='transcript_id').set_index('transcript_id')

    df = df.set_index('feature_id') 
    obs = df.columns
    obs = pd.DataFrame(obs, columns=['id']).set_index('id')


    adata = anndata.AnnData(obs=obs,
                            var=var,
                            X=df.T.values)
    
    adata.layers['total'] = adata.X
    adata.layers['total_raw'] = adata.X
    
    if save_path is not None:
        adata.write_h5ad(save_path)

    return(adata)

def filter_isoad(iso_ad, gene_ad, gene_name='gene_name',
                obs_copy=['leiden', 'state'],
                uns_copy=['leiden_colors', 'state_colors', 'neighbors', 'pca', 'umap'],
                obsp_copy=['connectivities', 'distances'],
                obsm_copy=['X_pca', 'X_umap'],
                save_path=None):
    iso_ad_ = iso_ad.copy()
    iso_ad_ = iso_ad_[np.intersect1d(iso_ad_.obs_names, gene_ad.obs_names)]
    gene_ad_ = gene_ad[iso_ad_.obs_names]
    iso_ad_ = iso_ad_[:, iso_ad.var[gene_name].isin(gene_ad.var[gene_name])]
    
    for obs in obs_copy:
        iso_ad_.obs[obs] = gene_ad_.obs[obs]
    for uns in uns_copy:
        iso_ad_.uns[uns] = gene_ad_.uns[uns]
    for obsp in obsp_copy:
        iso_ad_.obsp[obsp] = gene_ad_.obsp[obsp]
    for obsm in obsm_copy:
        iso_ad_.obsm[obsm] = gene_ad_.obsm[obsm]

    if save_path is not None:
        adata.write_h5ad(save_path)
    return iso_ad_



def to_rows_gff3(anno): #from pyranges
    rowdicts = []

    for line in list(anno):
        # stripping last white char if present
        lx = (it.split("=") for it in line.rstrip("; ").split(";"))
        rowdicts.append({k: v for k, v in lx})

    return pd.DataFrame.from_records(rowdicts).set_index(anno.index)

def get_interbed_df(intersect_bed, save_path=None):
    interbed = pd.read_csv(intersect_bed, sep='\t', names=range(15))
    extra = to_rows_gff3(interbed[14])
    interbed = pd.concat([interbed, extra], axis=1)
    interbed = interbed.drop(14, axis=1)
    interbed[1] += 1 #re-make bedfile into 1-based
    interbed.insert(0, 'exon_id_isoq',
                     [f"{l[0]}:{l[1]}-{l[2]}:{l[3]}:{l[5]}" for _, l in interbed.iterrows()]) #make index
    if save_path is not None:
        interbed.to_csv(save_path, index=None)

    return interbed

def filter_interbed_geneid(interbed):
    keep_index = []
    for _, line in interbed.iterrows():
        geneidx = line[5].split(',')
        geneidy = line['gene_id']
        for g in geneidx:
            if g == geneidy:
                keep_index.append(_)
                continue
    interbed_ = interbed.loc[keep_index]
    return interbed_

def get_overlap(interbed):
    interbed = interbed.copy()
    interbed['length_x'] = interbed[2] - interbed[1] +1  #one-based
    interbed['length_y'] = interbed[10] - interbed[9] +1
    
    overlap = []
    for _, line in interbed.iterrows():
        positions = [line[1], line[2], line[9], line[10]]
        positions_u = np.unique(positions)
        
        len_x = line['length_x']
        len_y = line['length_y']

        if line[1] == line[9]:
            if line[2] == line[10]:
                overlap.append(len_x)
            else:
                overlap.append(min(len_x, len_y))
        elif line[2] == line[10]:
            overlap.append(min(len_x, len_y))
        elif len(positions_u)==3:
            overlap.append(1)
        else:
            try:
                overlap.append(positions_u[2] - positions_u[1]+1)
            except:
                overlap.append(0)
    interbed['overlap'] = overlap
    interbed['overlap_prop_x'] = interbed['overlap']*100/interbed['length_x']
    interbed['overlap_prop_y'] = interbed['overlap']*100/interbed['length_y']
    return interbed

def filter_interbed_overlap(interbed, overlap_prop_x=100, overlap_prop_y=100):
    interbed_ = interbed[
    (interbed['overlap_prop_x']>=overlap_prop_x)&(interbed['overlap_prop_y']>=overlap_prop_y)
    ]
    return interbed_     


def get_adata_from_exon_isoquant(datapath, save_path=None):
    df = pd.read_csv(datapath, sep='\t')
    df_ = df.dropna(subset='group_id')

    df_include = df_.pivot_table(index=['chr', 'start', 'end', 'strand', 'flags', 'gene_ids'],
                             columns='group_id',
                             values='include_counts')
    df_exclude = df_.pivot_table(index=['chr', 'start', 'end', 'strand', 'flags', 'gene_ids'],
                             columns='group_id',
                             values='exclude_counts')
    def tream_pivot(df_c):
        df_c = df_c.fillna(0).reset_index()
        df_c['exon_id_isoq'] = [
            f"{l['chr']}:{l['start']}-{l['end']}:{l['strand']}:{l['gene_ids']}" for _, l in df_c.iterrows()]
        var = df_c[['chr', 'start', 'end', 'strand', 'flags', 'gene_ids', 'exon_id_isoq']]
        var = var.set_index('exon_id_isoq')
        var['length'] = var['end'] - var['start'] + 1
        df_c = df_c.set_index('exon_id_isoq')
        df_c = df_c.drop(['chr', 'start', 'end', 'strand', 'flags', 'gene_ids'],
                         axis=1)
        return(df_c, var)
    df_include, var = tream_pivot(df_include)
    df_exclude, var = tream_pivot(df_exclude)
    obs = df_include.columns
    obs = pd.DataFrame(obs, columns=['group_id']).set_index('group_id')

    
    df_sum = df_include + df_exclude

    adata = anndata.AnnData(obs=obs,
                            var=var,
                            X=df_include.T.values)
    adata.layers['total_raw'] = df_sum.T.values
    adata.layers['include_raw'] = df_include.T.values
    adata.layers['exclude_raw'] = df_exclude.T.values

    if save_path is not None:
        adata.write_h5ad(save_path)

    return(adata)

def filter_exonad(exon_ad, gene_ad,
                obs_copy=['leiden', 'state'],
                uns_copy=['leiden_colors', 'state_colors', 'neighbors', 'pca', 'umap'],
                obsp_copy=['connectivities', 'distances'],
                obsm_copy=['X_pca', 'X_umap'],
                save_path=None):
    exon_ad_ = exon_ad.copy()
    exon_ad_ = exon_ad_[np.intersect1d(exon_ad_.obs_names, gene_ad.obs_names)]
    gene_ad_ = gene_ad[exon_ad_.obs_names]
    
    for obs in obs_copy:
        exon_ad_.obs[obs] = gene_ad_.obs[obs]
    for uns in uns_copy:
        exon_ad_.uns[uns] = gene_ad_.uns[uns]
    for obsp in obsp_copy:
        exon_ad_.obsp[obsp] = gene_ad_.obsp[obsp]
    for obsm in obsm_copy:
        exon_ad_.obsm[obsm] = gene_ad_.obsm[obsm]

    if save_path is not None:
        adata.write_h5ad(save_path)
    return exon_ad_

    
    
    
    




