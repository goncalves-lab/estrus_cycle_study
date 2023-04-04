import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
import pandas as pd
from scvi.model import CondSCVI, DestVI
import anndata as ad
import random
import sys
import destvi_utils

sc.set_figure_params(figsize=(4, 4), frameon=False)



os.environ['HTTP_PROXY']="http://www-int.dkfz-heidelberg.de:80"
os.environ['HTTPS_PROXY']="http://www-int.dkfz-heidelberg.de:80"

if __name__ == '__main__':
    os.chdir(sys.argv[1])
    metadata = pd.read_csv("".join(["metadata_", sys.argv[2],".csv"]), sep=',', index_col=0)
    reads = pd.read_csv("".join(["expression_", sys.argv[2],".csv"]), sep=',',index_col=0)
    sc_adata=ad.AnnData(X=reads, obs=metadata)
    G = 2000
    sc.pp.filter_genes(sc_adata, min_counts=10)

    sc_adata.layers["counts"] = sc_adata.X.copy()

    sc.pp.highly_variable_genes(
        sc_adata,
        n_top_genes=G,
        subset=True,
        layer="counts",
        flavor="seurat_v3"
    )
    G = 2000
    sc.pp.filter_genes(sc_adata, min_counts=10)

    sc_adata.layers["counts"] = sc_adata.X.copy()

    sc.pp.highly_variable_genes(
        sc_adata,
        n_top_genes=G,
        subset=True,
        layer="counts",
        flavor="seurat_v3"
    )
    sp_metadata = pd.read_csv("".join(["metadata_sp_", sys.argv[3],".csv"]), sep=',', index_col=0)
    sp_reads = pd.read_csv("".join(["expression_sp_", sys.argv[3],".csv"]), sep=',',index_col=0)
    os.chdir(sys.argv[4])
    data = np.empty((len(sp_metadata), 2))
    data[:,0] = sp_metadata['imagerow']
    data[:,1] = sp_metadata['imagecol']

    st_adata=ad.AnnData(X=sp_reads, obs=sp_metadata)
    st_adata.layers["counts"] = st_adata.X.copy()
    st_adata.obsm['spatial'] = data
    sc.pp.normalize_total(st_adata, target_sum=10e4)
    sc.pp.log1p(st_adata)
    st_adata.raw = st_adata
    # filter genes to be the same on the spatial data
    intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
    st_adata = st_adata[:, intersect].copy()
    sc_adata = sc_adata[:, intersect].copy()
    G = len(intersect)
    CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key="cell_types")
    sc_model = CondSCVI(sc_adata, weight_obs=False)
    sc_model.train(max_epochs=500)
    f = sc_model.history["elbo_train"].iloc[5:]
    f.plot().get_figure().savefig("".join(["sc_training_", sys.argv[333],".pdf"]))
    DestVI.setup_anndata(st_adata, layer="counts")
    st_model = DestVI.from_rna_model(st_adata, sc_model, vamp_prior_p=5)
    st_model.train(max_epochs=5000)
    st_model.history["elbo_train"].iloc[10:].plot()
    f = st_model.history["elbo_train"].iloc[5:]
    f.plot().get_figure().savefig("".join(["st_training_", sys.argv[3],".pdf"]))
    st_adata.obsm["proportions"] = st_model.get_proportions()
    ct_list = ["APC", "BC", "EC","GC","SC","EpC","LC","TC"]
    for ct in ct_list:
        data = st_adata.obsm["proportions"][ct].values
        st_adata.obs[ct] = np.clip(data, 0, np.quantile(data, 0.99))
    sc.pl.embedding(st_adata, basis="spatial", color=ct_list, cmap='Reds', s=80, save="".join(["proportions_",sys.argv[3],".pdf"]))
    st_adata.write("".join(["st_adata_",sys.argv[3],".h5ad"]))
    
    
    
