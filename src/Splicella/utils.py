import pandas as pd
import scanpy as sc
import numpy as np
import pyranges as pr
import anndata

def get_annotr_fromgtf(gtf_path, save_path=None):
    anno = pr.read_gtf('refs/gencode.v48.primary_assembly.annotation.gtf.gz').df
    anno_tr = anno[anno['Feature'] == 'transcript']
    anno_tr_ = anno_tr[['transcript_id', 'gene_name', 'gene_id',
                    'gene_type', 'transcript_type']]
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


def get_adata_from_exon(filepath, save_path):
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
        var['length'] = var['end'] - var['start'] + 1 #one-based
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
    

def filter_exonad(exon_ad, gene_ad, gene_name='gene_name',
                obs_copy=['leiden', 'state'],
                uns_copy=['leiden_colors', 'state_colors', 'neighbors', 'pca', 'umap'],
                obsp_copy=['connectivities', 'distances'],
                obsm_copy=['X_pca', 'X_umap'],
                save_path=None):
    iso_ad_ = iso_ad.copy()
    iso_ad_ = iso_ad_[gene_ad.obs_names]

    for obs in obs_copy:
        iso_ad_.obs[obs] = gene_ad.obs[obs]
    for obsp in obsp_copy:
        iso_ad_.obs[obsp] = gene_ad.obsp[obsp]
    for obsm in obsm_copy:
        iso_ad_.obs[obsm] = gene_ad.obsp[obsm]

    if save_path is not None:
        adata.write_h5ad(save_path)
    return iso_ad_ 
    

def get_adata_from_exon(filepath, overlap_filter=True, overlap_x=100, overlap_y=100, save_path=None):
    def to_rows_gff3(anno):
        rowdicts = []
    
        for line in list(anno):
            # stripping last white char if present
            lx = (it.split("=") for it in line.rstrip("; ").split(";"))
            rowdicts.append({k: v for k, v in lx})
    
        return pd.DataFrame.from_records(rowdicts).set_index(anno.index)
    
    df = pd.read_csv('data/RCC1/scNanoGPS/exon_intersect_wab.bed', sep='\t',
                      names=range(15), )
    df[1] +=1 #as scNanoGPS output is 1-based and bed file is 0 based
    
    #make attributes as multiple column
    extra = to_rows_gff3(df[14])
    df = pd.concat([interbed, extra], axis=1)

    
    keep_index = []
    for _, line in df.iterrows():
        geneidx = line[5].split(',')
        geneidy = line['gene_id']
        for g in geneidx:
            if g == geneidy:
                keep_index.append(_)
                continue

    df_ = df.loc[keep_index]
    df_['length_x'] = df_[2] - df_[1] +1
    df_['length_y'] = df_[10] - df_[9] +1

    overlap = []
    for _, line in df_.iterrows():
        positions = [line[1], line[2], line[9], line[10]]
        positions_u = np.unique(positions)
        if len(positions_u)==2:
            overlap.append(positions_u[1] - positions_u[0])
        elif len(positions_u)==3:
            if positions[0] == positions[2]:
                overlap.append(positions_u[1] - positions_u[0])
            else:
                overlap.append(positions_u[2] - positions_u[1])
        else:
            overlap.append(positions_u[2] - positions_u[1])

    df_['overlap'] = overlap
    df_['overlap_prop_x'] = df_['overlap']*100/df_['length_x']
    df_['overlap_prop_y'] = df_['overlap']*100/df_['length_y']

    df_.insert(0, 'exon_id_isoq',
                        [f"{l[0]}:{l[1]}-{l[2]}:{l[3]}:{l[5]}" for _, l in interbed_.iterrows()]) 

    if overlap_filter:
        keeploc = []
        for eid, df in interbed_.groupby('exon_id_isoq'):
            if len(df[(df['overlap_prop_x']==100)&(df['overlap_prop_y']==100)])>0:
                df = df[(df['overlap_prop_x']==100)&(df['overlap_prop_y']==100)]
                keeploc.extend(df.index.to_list())
            else:
                continue
        df_ = df_.loc[keeploc].sort_index()
        
    if save_path is not None:
        df_.to_csv(save_path)
        
    return df_
    
    




