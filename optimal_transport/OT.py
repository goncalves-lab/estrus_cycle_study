#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#############
#this script uses files generated in the script nmf_preparation.R and nmf matrix generated using PLANC library 

import numpy as np
import scipy
from sklearn.metrics import pairwise_distances
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import pdist
import os
import pandas as pd
import anndata as ad
import sys
from scipy.sparse import vstack


def group_distance(X: np.ndarray,
                   idx0: np.ndarray,
                   idx1: np.ndarray,
                   metric: str = 'munkres_cost',
                   n_jobs: int = 1) -> float:
    '''Calculate the distance between two groups of points in `X`
    Parameters
    ----------
    X : np.ndarray
        [Cells, Features] embedding.
    idx0 : np.ndarray
        [Cells,] index of cells in group 0.
    idx1 : np.ndarray
        [Cells,] index of cells in group 1.
    metric : str
        distance metric to use.
        'centroid_distance_euclidean' - euclidean dist between centroids.
        'centroid_distance_manhattan' - manhattan dist between centroids.
        'centroid_distance_cosine' - cos dist between centroids.
        'average' - mean pairwise euclidean dist between members of 
            group `0` and members of group `1`.
        'median' - median pairwise euclidean dist between members of 
            group `0` and members of group `1`.
        'cosine' - mean pairswise cosine dist between members of 
            group `0` and members of group `1`.
        'munkres_cost' - munkres min cost solution for optimal transport
            of group `0` to group `1`.
    Returns
    -------
    dist : float
        distance between groups 0 and 1.
    '''
    # dm[i, j] is distance between X[i,:] and X[j,:]

    if metric.lower() == 'munkres_cost':

        x_0 = X[idx0, :]
        x_1 = X[idx1, :]
        dm = pairwise_distances(x_0, x_1, n_jobs=n_jobs, metric='euclidean')
        row_idx, col_idx = linear_sum_assignment(dm)
        dist = np.sum(dm[row_idx, col_idx])

    else:
        raise ValueError('arg to `metric` %s not understood' % metric)

    return dist

def random_sample(groupbys: list,
                  groups: list,
                  n: int = 200) -> np.ndarray:
    '''Generate a random index of cells where cells 
    are drawn from `groups` in `groupbys`. Multiple `groupbys`
    can be provided simultaneously for expedient random
    sampling of subpopulations.
    Parameters
    ----------
    groupbys : list
        list of [Cells,] vectors of group assignments for 
        each group to sample from.
        e.g. [
            np.array(['typeA', 'typeA', 'typeB']),
            np.array(['0', '1', '0']),
        ]
        where each of the arrays is a unique grouping variable.
    groups : list
        list of str values in groupbys to sample from.
        e.g. ['typeA', '0'] would sample only the first
        cell from the above `groupbys` example.
    n : int
        sample size.
    Returns
    -------
    idx : np.ndarray
        [n,] index of integers.
    Notes
    -----
    If sample size `n` is larger than the number of cells in a group,
    we reduce the sample size to `0.8*n_cells_in_group`.
    Examples
    --------
    >>> groupbys = [
    ...   adata.obs['cell_type'],
    ...   adata.obs['louvain'],
    ... ]
    >>> # sample only B cells in cluster '3'
    >>> groups = ['B cell', '3']
    >>> # returns indices of random B cells
    >>> # in louvain cluster '3'
    >>> idx = random_sample(groupbys, groups)
    '''
    assert len(groupbys) == len(groups)
    bidxs = []
    for i in range(len(groups)):
        group_bidx = (groupbys[i] == groups[i])
        bidxs.append(group_bidx)
    bidx = np.logical_and.reduce(bidxs)
    group_idx = np.where(bidx)[0].astype(np.int)

    if int(0.8*len(group_idx)) < n:
        #print(f'Cannot draw {n} cells from group of size {len(group_idx)}.')
        n = 0.8*len(group_idx)
        #print(f'Reducing sample size to {n}.')

    idx = np.random.choice(group_idx,
                           size=int(n),
                           replace=False).astype(np.int)
    return idx


#calculating OT distance for all possible combinations estrus phases for all cell types

def distance_calculation_final(days, cells, adata, num_set):
    row_list=[]
    #final_list=[]
    for cell in cells:
        for day in days:
            dist_0_1_average=[]
            print(cell)
            print(day)
            for counter in range(100):
                idx_0 = random_sample(
                groupbys = [
                    adata.obs['day'], 
                    adata.obs['cells']
                    ],
                    groups = [day[0], cell],
                    n=num_set,
                    )
                idx_1 = random_sample(
                    groupbys = [
                    adata.obs['day'], 
                    adata.obs['cells']
                    ], 
                    groups = [day[1], cell],
                    n=num_set,
                    )
                
                   #print(len(idx_0), len(idx_1))
                if len(idx_0) > 0 and len(idx_1) > 0:
                      dist_0_1 = group_distance(
                      X=adata.X,
                                idx0=idx_0,
                                idx1=idx_1,
                                )
                      dist_0_1_average.append(dist_0_1)
                            #print(dist_0_1_average)
            if len(dist_0_1_average) > 0:
                dist_0_1_final = sum(dist_0_1_average) / len(dist_0_1_average)
                name=[cell, str(day[0]), str(day[1])]
                name=name="_".join(name)
                dict_dist= {
                                    "condition": name, 
                                    "distance": dist_0_1_final,
                                    }
                row_list.append(dict_dist)
                
                #dict_dist={
                #    "conditon":row_list[0]["condition"],
                #    "distance":sum(elem['distance'] for elem in row_list) / len(row_list)
                #    }
                #final_list.append(dict_dist)
            
        #print(row_list)
    df_distance=pd.DataFrame(row_list)
    return df_distance



#loading all necessary files
metadata = pd.read_csv("metadata_hvg_agg", sep=',', index_col=0)
reads = pd.read_csv("nmf_estrus", sep=',',index_col=0)
reads=scipy.sparse.csr_matrix(reads)

adata_estrus=ad.AnnData(X=reads, obs=metadata)

cells=adata_estrus.obs['cells'].to_list()
cells=set(cells)
cells=list(cells)

days=[("proestrus", "estrus"), ("estrus", "metestrus") ,("metestrus", "diestrus"), ("diestrus", "proestrus")]

df_distance=distance_calculation_final(days, cells, adata_estrus, 100)

df_distance.to_csv("df_distance",index=False)




     










