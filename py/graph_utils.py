import os 
import pandas as pd
import numpy as np
import networkx as nx
import itertools

# Utility library to deal with graphs, neighbors and similar

def get_children_w_graph(treats, total_treats,G, Zs, data_cols):

    # select unique neighbors
    adj = [[n for n in G.neighbors(t.split("_")[0])] for t in treats if G.has_node(t.split("_")[0])]
    adj = [item for sublist in adj for item in sublist]
    adj.sort()
    adj = list(k for k,_ in itertools.groupby(adj))

    # add all their "variants" (somatic, LOH, ...)
    adj_ext = [[n for n in total_treats if "_"+a+"_" in "_"+n] for a in adj]
    adj_ext = [item for sublist in adj_ext for item in sublist]


    # cleanup
    rem_list = set([a for a in adj_ext if a in treats or a in Zs or a not in data_cols])
    adj_ext = [a for a in adj_ext if a not in rem_list]

    return adj_ext


def calculate_children_number_w_graph(total_treats,G,Zs,df, lv_stop_calculus):
    res = [len(total_treats)]
    prev = total_treats

    for lv in range(1,lv_stop_calculus+1): # lv counts the number of Xs

        curr_lev = []
        for treats in prev:
            if lv == 1:
                treats = [treats]

            adj = get_children_w_graph(treats, total_treats,G,Zs, data_cols = df.columns.tolist())

            for v in adj:
                child_treats = sorted( treats + [v])
                curr_lev.append(child_treats)
        
        rem_dup = []
        
        curr_lev.sort()
        rem_dup = list(k for k,_ in itertools.groupby(curr_lev))

        res.append(len(rem_dup))
        prev = rem_dup

    return res


