import os 
import pandas as pd
import numpy as np
from scipy.spatial import distance

# Library used for the cleanup proposed after finding the best rule, for discovering rules that effectively map all the rule space

def remove_genes_similar_to_solution(best_node_state, df, total_treats, thr = 0.01):
    data_path, cfnds, treats, treats_vals, outcome, outcome_vals, per_strata_data = best_node_state
    rem_list = []

    m = [[distance.cityblock(df[x],df[y])/len(df) for x in total_treats] for y in treats]

    for i in range(len(treats)):
        rem_list = rem_list + [total_treats[j] for j in range(len(total_treats)) if m[i][j]<thr or m[i][j]>1-thr]
        # print([m[i][j] for j in range(len(total_treats)) if m[i][j]<thr or m[i][j]>1-thr])

    return [x for x in total_treats if x not in rem_list]

