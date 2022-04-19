#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 19:06:52 2021


@author: molbio
"""
# to run this script following files are necessary: list of TFs, motifs and cis_Target_db (obtained from SCENIC depository), and 
# normalized read matrix
import os
import glob
import sys
import socket
import subprocess
import pickle

import pandas as pd
import numpy as np
#import scipy
#import seaborn as sns

#import anndata as ad
#import loompy as lp

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, save_to_yaml, load_from_yaml
from pyscenic import prune
from pyscenic.prune import prune2df
from pyscenic.transform import df2regulons      #was missing - is this right?
from pyscenic.aucell import aucell
from pyscenic.rss import regulon_specificity_scores

from dask.distributed import Client, LocalCluster   #needed to run grnboost2


if __name__ == '__main__':
    HOSTNAME = socket.gethostname()
    DKFZ_CLUSTER = HOSTNAME != 'molbio'
    if DKFZ_CLUSTER:
        DATA_ROOT = '/internal/'
        DATA_BRIDGE = 'scenic'
        DATA_TAG = '/estrus_cycle/'
    else:
        DATA_ROOT = os.path.expanduser('~/')
        DATA_BRIDGE = 'scenic_cluster'
        DATA_TAG = ''
    wdir = f"{DATA_ROOT}{DATA_TAG}/{DATA_BRIDGE}"
    os.chdir(wdir)
    f_tfs = f"{wdir}/allTFs_mm.txt"
    reads = pd.read_csv(sys.argv[1], sep=',',index_col=0)
    f_db_glob =  f"{wdir}/cis_Target_db/*feather"
    f_motif_path = f"{wdir}/cis_Target_db/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
        # motif databases (e.g. "cis_Target_db/motifs-v9-nr.mgi-m0.001-o0.0.tbl")
        #file contained in wd
        #load ranking databases
    db_fnames = glob.glob(f_db_glob)
    def name(fname):
        return os.path.basename(fname).split(".")[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
    tf_names = load_tf_names(f_tfs)
    local_directory = f'{os.getenv("SCRATCHDIR")}/{os.getenv("LSB_JOBID")}'
    local_cluster = LocalCluster(
        n_workers=20,
    local_directory=local_directory)
    client = Client(local_cluster)
    adjacencies = grnboost2(expression_data=reads,
                                tf_names=tf_names, seed = 123, verbose=True, client_or_address=client)
    adjacencies.to_csv(sys.argv[2], sep = '\t')
    #adjacencies  = pd.read_csv(sys.argv[2], sep='\t', header=0, index_col=0)
    modules = list(modules_from_adjacencies(adjacencies, reads))
    df = prune2df(dbs, modules, f_motif_path, client_or_address="dask_multiprocessing") #originally "local"
    regulons = df2regulons(df)
    regulons_df = pd.DataFrame({
        'tf': [regulon.name for regulon in regulons],
        'targets': [
            ' '.join(list(regulon.gene2weight)) 
            for regulon in regulons
        ]
    })
    regulons_df.to_csv(sys.argv[3], index=False)
    auc_mtx = aucell(reads, regulons, num_workers=20)
    auc_mtx.to_csv(sys.argv[4], sep = '\t')
