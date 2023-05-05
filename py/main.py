from re import X 
import numpy as np
import pandas as pd
from typing import Dict, List
from allstar import BnB_effect_calculator, BnB_effect_calculator_naive
import pybnb
import os
import networkx as nx
import graph_utils
from math import comb
import argparse
import ds_cleanup_utils


parser=argparse.ArgumentParser('main', description="Launch ALLSTAR, the corrected and reliable causal rule discovery tool for high-throughput datasets", 
                                       epilog='Read more at the link: *WILL BE ADDED AFTER REVISION*')
parser.add_argument('-f', '--data_path', required=True, help='Path of the dataset to be investigated. It must follow the following block order: confounders Z, treatments X, output Y')
parser.add_argument('-s', '--split_column', type=int, required=True, help='Number of confounders. It is used internally to subdivide the dataset into the founding blocks')
parser.add_argument('-w', '--wall_time_limit', type=int, default = 60*10, help='Time limit to reach convergence for each calculation. At the end of this time the best rule found will be accounted')
parser.add_argument('-v', '--version', type=int, default = 3, required=True, help="1 = ")
parser.add_argument('-t', '--thr', type=int, default=0.01, help='Threshold of similarity between treatments. After the best rule is found, ALLSTAR removes treatments that do not differ more than this threshold from the ones found within the rule')
parser.add_argument('--max_rule_len', type=int, default=4, help='Maximum number of confounders + treatments to be accounted when searching for the best rule')
parser.add_argument('--graph_path', default = 'PPI/FIsInGene_122921_with_annotations.txt', help='Path of the interaction graph between treatments')
parser.add_argument('--alpha_set', type=float, default=0.05, help='Confidence level of the algorithm')
parser.add_argument('--use_MHT_correction', action='store_true', help='If selected, ALLSTAR corrects for Multiple Hypothesis Testing')
parser.add_argument('--use_graph_per_children', action='store_true', help='If selected, ALLSTAR restricts the search to rules including treatments strictly connected within the interaction graph')
parser.add_argument('--include_unmutated_treats', action='store_true', help='If selected, ALLSTAR includes the untreated values in the rule discovery process')

args=parser.parse_args()

data_path = args.data_path
split_column = args.split_column
wall_time_limit = args.wall_time_limit
max_rule_len = args.max_rule_len
alpha_set = args.alpha_set
use_MHT_correction = args.use_MHT_correction
use_graph_per_children = args.use_graph_per_children
include_unmutated_treats = args.include_unmutated_treats
graph_path = args.graph_path
version = args.version
thr = args.thr

mod_version = version
if version > 100:
    Zs = []
    mod_version = version - 100

#"results/dataset_name/v2"
if not os.path.exists("results"):
    os.mkdir("results")
if mod_version < 10:
    dataset_name_only = data_path.split("/")[-1].split(".")[0]
else:
    dataset_name_only = "synth_v"+str(mod_version)
if not os.path.exists("results/"+dataset_name_only):
    os.mkdir("results/"+dataset_name_only)
if not os.path.exists("results/ALLSTAR_runs"):
    os.mkdir("results/ALLSTAR_runs")
if not os.path.exists("results/"+dataset_name_only+"/v"+str(mod_version)):
    # shutil.rmtree("results/"+dataset_name_only+"/v"+str(mod_version))
    os.mkdir("results/"+dataset_name_only+"/v"+str(mod_version))
output_dir = "results/"+dataset_name_only+"/v"+str(mod_version)

def get_file_name(output_dir,y_vals,cfnd,include_unmutated_treats,use_MHT_correction,use_graph_per_children, prefix = ""):
    output_file_name = output_dir + "/" + prefix + "Y_in_" + str(y_vals)
    if len(cfnd)>0:
        output_file_name += "|"+str(cfnd).replace("[","").replace("]","")
    output_file_name += "_max_len_"+str(max_rule_len)
    if include_unmutated_treats:
        output_file_name += "_w_0s"
    if use_MHT_correction:
        output_file_name += "_corr"
    if use_graph_per_children:
        output_file_name += "_graph"
    output_file_name +=".txt"
    return output_file_name


if True: # default parameters setting
    data = pd.read_csv(data_path)
    data_colnames = data.columns.values.tolist()
    Y_variable = data_colnames[-1]
    # Drop last column as it is the outcome 
    data_colnames = data_colnames[:-1]
    # Isolate treatments
    Xs = data_colnames[split_column:]
    Zs = data_colnames[:split_column]

    # Calculate effects
    cfnds = []
    treats = []
    treats_vals = [[]]


    Y_classes = [[a] for a in sorted(data[Y_variable].unique())]
    Y_values = len(Y_classes)

if version > 100:
    Zs = []
    version = version - 100

### modifica valori input
if version == 1:  # default BnB normale + correzioni
    use_graph_per_children = False
if version == 2:  # BnB grafo + correzioni
    pass
if version == 3:  # BnB "ricorrente" che elimina le vars
    pass
if version == 11: # synthetic G1: dimensionality
    use_graph_per_children = False
    use_MHT_correction = False
if version == 12: # synthetic G2: graph utility
    use_MHT_correction = True
if version == 13: # synthetic G3: heterogeneous iterative search
    use_MHT_correction = True
    use_graph_per_children = True
    Y_classes = [[1]]
    Y_values = len(Y_classes)



assert np.min(data[Xs].values) >= 0 and np.max(data[Xs].values)<=1, "Treatments do not take values in 0,1 only"

if True: # calculus of correction for MHT
    G = None
    if use_graph_per_children:
        df_graph = pd.read_table(graph_path)
        G = nx.from_pandas_edgelist(df_graph, 'Gene1', 'Gene2', 'Score')
    ########### NB per i figli inserisco solo i valori in Xs.
    ########### se voglio quindi tenere o togliere qualcosa, conviene farlo PRIMA di dividere la matrice
    ########### nelle varie X Y e Z
    if use_MHT_correction:
        print("- Starting MHT correction calculus")
        lxs = len(Xs)
        if use_graph_per_children:
            number_hyp_per_level = graph_utils.calculate_children_number_w_graph(Xs,G,Zs,data, max_rule_len-1)  # non cumulative number of distinct hypotheses per level if using the Protein-Protein graph
        else:
            # if include_unmutated_treats:
            #     number_hyp_per_level = [(2**j)*comb(lxs,j) for j in range(1,max_rule_len)]
            # else:            
            number_hyp_per_level = [comb(lxs,i) for i in range(1,max_rule_len)]  # non cumulative number of distinct hypotheses per level if using the Protein-Protein graph

        if include_unmutated_treats:
            number_hyp_per_level = [(2**j)*number_hyp_per_level[j-1] for j in range(1,max_rule_len)]

        stop_level = max_rule_len-1
        cumulative_number_of_Xs = sum(number_hyp_per_level[:stop_level])
        # print(cumulative_number_of_Xs)
        val_zs = sum([len(set(data[z_name].values)) for z_name in Zs])
        # print(Zs,split_column)
        # print(val_zs,[len(set(data[z_name].values)) for z_name in Zs])
        tot_N_hyps = cumulative_number_of_Xs * Y_values * (val_zs + 1)  
        # print(tot_N_hyps)

        alpha_set = alpha_set/tot_N_hyps
        print("- Ended MHT correction calculus")


### multiple runs on real world datasets

if version == 1:  # default BnB standard + MHT correction
    for cfnd in [[]]+Zs:
    # for cfnd in Zs:
        if len(cfnd)!=0:
            cfnd = [cfnd]
        for y_v in Y_classes:
            output_file_name = get_file_name(output_dir,y_v,cfnd,include_unmutated_treats,use_MHT_correction,use_graph_per_children)
            state = BnB_effect_calculator(data_path, cfnd, treats, treats_vals, Xs, Y_variable, y_v, 
                wall_time_limit, max_rule_len, alpha_set, use_graph_per_children, include_unmutated_treats, G, output_file_name)

if version == 2:  # BnB + graph + MHT correction
    for cfnd in [[]]+Zs:
        if len(cfnd)!=0:
            cfnd = [cfnd]
        for y_v in Y_classes:
            output_file_name = get_file_name(output_dir,y_v,cfnd,include_unmutated_treats,use_MHT_correction,use_graph_per_children)
            state = BnB_effect_calculator(data_path, cfnd, treats, treats_vals, Xs, Y_variable, y_v, 
                wall_time_limit, max_rule_len, alpha_set, use_graph_per_children, include_unmutated_treats, G, output_file_name)

if version == 3:  # BnB with iterative cleaning
    n_iter = 3
    for cfnd in [[]]+Zs:
        if len(cfnd)!=0:
            cfnd = [cfnd]
        for y_v in Y_classes:
            Xs = data_colnames[split_column:]
            output_file_name = get_file_name(output_dir,y_v,cfnd,include_unmutated_treats,use_MHT_correction,use_graph_per_children)
            for i in range(n_iter):
                    state = BnB_effect_calculator(data_path, cfnd, treats, treats_vals, Xs, Y_variable, y_v, 
                        wall_time_limit, max_rule_len, alpha_set, use_graph_per_children, include_unmutated_treats, G, output_file_name)

                    # update treats
                    Xs = ds_cleanup_utils.remove_genes_similar_to_solution(state, data, Xs, thr = thr)

### multiple runs on synthetic datasets
if version == 11:  # synthetic G1: dimensionality 

    dataset_name_only = data_path.split("/")[-1].split(".")[0]

    # naive
    for cfnd in [[]]+Zs:
        if len(cfnd)!=0:
            cfnd = [cfnd]
        for y_v in Y_classes:
            output_file_name = get_file_name(output_dir,y_v,cfnd,include_unmutated_treats,use_MHT_correction,use_graph_per_children,prefix = "G1_naive_"+dataset_name_only)
            state = BnB_effect_calculator_naive(data_path, cfnd, treats, treats_vals, Xs, Y_variable, y_v, 
                wall_time_limit, max_rule_len, include_unmutated_treats, output_file_name)
    # reliable
    for cfnd in [[]]+Zs:
        if len(cfnd)!=0:
            cfnd = [cfnd]
        for y_v in Y_classes:
            output_file_name = get_file_name(output_dir,y_v,cfnd,include_unmutated_treats,use_MHT_correction,use_graph_per_children,prefix = "G1_reliable_"+dataset_name_only)
            state = BnB_effect_calculator(data_path, cfnd, treats, treats_vals, Xs, Y_variable, y_v, 
                wall_time_limit, max_rule_len, alpha_set, use_graph_per_children, include_unmutated_treats, G, output_file_name)
    # corrected
    use_MHT_correction = True
    
    print("- Starting MHT correction calculus")
    lxs = len(Xs)
    number_hyp_per_level = [comb(lxs,i) for i in range(1,max_rule_len)]  # non cumulative number of distinct hypotheses per level if using the Protein-Protein graph
    stop_level = max_rule_len-1
    cumulative_number_of_Xs = sum(number_hyp_per_level[:stop_level])

    tot_N_hyps = cumulative_number_of_Xs * Y_values * (len(Zs) + 1)  
    
    if include_unmutated_treats:
        tot_N_hyps = tot_N_hyps * 2

    alpha_set = alpha_set/tot_N_hyps
    print("- Ended MHT correction calculus")

    for cfnd in [[]]+Zs:
        if len(cfnd)!=0:
            cfnd = [cfnd]
        for y_v in Y_classes:
            output_file_name = get_file_name(output_dir,y_v,cfnd,include_unmutated_treats,use_MHT_correction,use_graph_per_children,prefix = "G1_corrected_"+dataset_name_only)
            state = BnB_effect_calculator(data_path, cfnd, treats, treats_vals, Xs, Y_variable, y_v, 
                wall_time_limit, max_rule_len, alpha_set, use_graph_per_children, include_unmutated_treats, G, output_file_name)

if version == 12:  # synthetic G2: graph utility 

    dataset_name_only = data_path.split("/")[-1].split(".")[0]

    for cfnd in [[]]+Zs:
        if len(cfnd)!=0:
            cfnd = [cfnd]
        for y_v in Y_classes:
            output_file_name = get_file_name(output_dir,y_v,cfnd,include_unmutated_treats,use_MHT_correction,use_graph_per_children,prefix = "G2_"+dataset_name_only)
            state = BnB_effect_calculator(data_path, cfnd, treats, treats_vals, Xs, Y_variable, y_v, 
                wall_time_limit, max_rule_len, alpha_set, use_graph_per_children, include_unmutated_treats, G, output_file_name)

if version == 13:  # synthetic G3: heterogeneous iterative search 

    dataset_name_only = data_path.split("/")[-1].split(".")[0]
    keep_Xs = Xs

    n_iter = 4
    for cfnd in [[]]+Zs:
        if len(cfnd)!=0:
            cfnd = [cfnd]
        for y_v in Y_classes:
            Xs = data_colnames[split_column:]
            output_file_name = get_file_name(output_dir,y_v,cfnd,include_unmutated_treats,use_MHT_correction,use_graph_per_children,prefix = "G3_light_"+dataset_name_only)
            for i in range(n_iter):
                    state = BnB_effect_calculator(data_path, cfnd, treats, treats_vals, Xs, Y_variable, y_v, 
                        wall_time_limit, max_rule_len, alpha_set, use_graph_per_children, include_unmutated_treats, G, output_file_name)

                    # update treats
                    Xs = ds_cleanup_utils.remove_genes_similar_to_solution(state,data,Xs, thr = 0.001)

    Xs = keep_Xs
    n_iter = 4
    for cfnd in [[]]+Zs:
        if len(cfnd)!=0:
            cfnd = [cfnd]
        for y_v in Y_classes:
            Xs = data_colnames[split_column:]
            output_file_name = get_file_name(output_dir,y_v,cfnd,include_unmutated_treats,use_MHT_correction,use_graph_per_children,prefix = "G3_strong_"+dataset_name_only)
            for i in range(n_iter):
                    state = BnB_effect_calculator(data_path, cfnd, treats, treats_vals, Xs, Y_variable, y_v, 
                        wall_time_limit, max_rule_len, alpha_set, use_graph_per_children, include_unmutated_treats, G, output_file_name)

                    # update treats
                    Xs = ds_cleanup_utils.remove_genes_similar_to_solution(state,data,Xs, thr = 0.05)
