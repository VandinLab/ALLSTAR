from re import X 
import sys
import pandas as pd
from typing import Dict, List
from effect_calculator import EffectCalculator
import bnb_methods
import pybnb
import networkx as nx
import graph_utils
from math import comb
import ds_cleanup_utils
import argparse


def BnB_effect_calculator(data_path, cfnds, treats, treats_vals, total_treats, Y_variable, Y_value, 
    wall_time_limit, max_rule_len, alpha_set, use_graph_per_children, include_unmutated_treats, G, output_file_name):

    try:
        orig_stdout = sys.stdout
        with open(output_file_name, 'a') as sys.stdout:
            problem = bnb_methods.ReliablePaperBnB(data_path, cfnds, treats, treats_vals, Y_variable, Y_value, total_treats, 
                max_rule_len = max_rule_len, alpha = alpha_set, use_graph_per_children = use_graph_per_children, 
                include_unmutated_treats = include_unmutated_treats, G=G)
            
            solver = pybnb.Solver()
            results = solver.solve(problem,
                absolute_gap = 1e-5,
                queue_strategy = "breadth",
                time_limit = wall_time_limit
                )

            if solver.is_dispatcher:            
                print(results.best_node.state) 

        sys.stdout = orig_stdout
        return results.best_node.state

    except OverflowError:
        print("Overflow Error TOT", sys.getsizeof(bnb_methods.Simple.SAVED_STATE))

def BnB_effect_calculator_naive(data_path, cfnds, treats, treats_vals, total_treats, Y_variable, Y_value, 
    wall_time_limit, max_rule_len, include_unmutated_treats, output_file_name):

    try:
        orig_stdout = sys.stdout
        with open(output_file_name, 'a') as sys.stdout:
            problem = bnb_methods.PaperBnB(data_path, cfnds, treats, treats_vals, Y_variable, Y_value, total_treats, 
                max_rule_len = max_rule_len, include_unmutated_treats = include_unmutated_treats)
            
            solver = pybnb.Solver()
            results = solver.solve(problem,
                absolute_gap = 1e-5,
                queue_strategy = "breadth",
                time_limit = wall_time_limit
                )

            if solver.is_dispatcher:            
                print(results.best_node.state) 

        sys.stdout = orig_stdout
        return results.best_node.state

    except OverflowError:
        print("Overflow Error TOT", sys.getsizeof(bnb_methods.Simple.SAVED_STATE))


def ALLSTAR_rule_calculus(data_path, X, Z, y_value, Y_variable, l, alpha, G, k, t, output_file_name):
    # inputs: alterations variables X, confounders variables Z, value y_value of target Y_variable , max. rule length l, confidence alpha, PPI graph G = (X, E), integer k, clean-up threshold t
    # output: at most k reliable rules with positive effect

    data = pd.read_csv(data_path)
    
    ### Calculating number N of hypotheses tested and correcting the confidence treshold aplha
    number_hyp_per_level = graph_utils.calculate_children_number_w_graph(X, G, Z, data, l-1)  # non cumulative number of distinct hypotheses per level if using the Protein-Protein graph
    number_hyp_per_level = [(2**j)*number_hyp_per_level[j-1] for j in range(1,l)] # Including also mutations equal 0 to the hypothesis count
    stop_level = l-1
    cumulative_number_of_Xs = sum(number_hyp_per_level[:stop_level])
    if Z != []:
        val_zs = sum([len(set(data[z_name].values)) for z_name in Z])  # need to correct for all possible Z values since I am marginalizing over Z, and therefore I infer all the different p(Y|sigma,Z=z) for all z values of Z
    else:
        val_zs = 1
    tot_N_hyps = cumulative_number_of_Xs * val_zs  
    alpha_set = alpha/tot_N_hyps
    
    # starting the (at most) k iteractions of branch-and-bound algorithm
    output = []
    max_rule_len = l + len(Z) # due to implementation choice, the length of a rule in the code is defined as the number of treatments + number of confounders, therefore, to have it alligned with the pseudocode, the max rule len is ell + size(Z)
    for i in range(k):
        # find the i-th rule
        state = BnB_effect_calculator(data_path, cfnds = Z, treats = [], treats_vals = [[]], total_treats = X, Y_variable = Y_variable, Y_value = y_value, wall_time_limit = None, max_rule_len = max_rule_len, alpha_set = alpha_set, use_graph_per_children = True, include_unmutated_treats = True, G = G, output_file_name = output_file_name)
        # update Xs
        X = ds_cleanup_utils.remove_genes_similar_to_solution(state, data, total_treats = X, thr = t)


def ALLSTAR_rule_parser(output_file_name):
    output = []
    with open(output_file_name, 'r') as f:
        data = f.readlines()
        for i in range(len(data)):
            if " - objective: " in data[i]:
                obj = float(data[i].replace(" - objective: ",""))
            if " - best_node: " in data[i]: # after 2 lines I have the rule encoding
                rule_enc = eval(data[i+2])
                rule = EffectCalculator.rule_tostring(cfnds = rule_enc[1], treats = rule_enc[2], treats_vals = rule_enc[3], outcome = rule_enc[4], outcome_vals = rule_enc[5])
                if obj > 0:
                    output.append((rule,obj))
    return output


parser=argparse.ArgumentParser('allstar', description="Launch ALLSTAR, the corrected and reliable causal rule discovery tool for high-throughput datasets", 
                                       epilog='Read more at the link: *WILL BE ADDED AFTER REVISION*')
parser.add_argument('-f', '--data_path', default='youforgotyour.csv', help='Path of the dataset to be investigated. It must follow the following block order: confounders Z, treatments X, output Y')
parser.add_argument('-s', '--split_column', type=int, default=99999, help='Number of confounders. It is used internally to subdivide the dataset into the founding blocks')
parser.add_argument('-k', '--top_k_rules', type=int, default=5, help='Number of rules to be returned.')
parser.add_argument('-t', '--thr', type=float, default=0.01, help='Threshold of similarity between treatments. After the best rule is found, ALLSTAR removes treatments that do not differ more than this threshold from the ones found within the rule')
parser.add_argument('-l', '--max_rule_len', type=int, default=4, help='Maximum number of confounders + treatments to be accounted when searching for the best rule')
parser.add_argument('-y', '--y_value', type=int, help='Value y to investigate. Choose one possible outcome of your target variable Y.')
parser.add_argument('-g', '--graph_path', default = 'PPI/FIsInGene_122921_with_annotations.txt', help='Path of the interaction graph between treatments')
parser.add_argument('-a', '--alpha', type=float, default=0.05, help='Confidence level of the algorithm')
parser.add_argument('--calculus_mode', action='store_true', help='If selected, ALLSTAR is set to work in causal inference modality.')
parser.add_argument('--parse_mode', action='store_true', help='If selected, ALLSTAR is set to parse results from a given result sheet.')
parser.add_argument('-o', '--output_file_name', default='output.txt', help='Output file to be parsed')

args=parser.parse_args()

data_path = args.data_path
split_column = args.split_column
max_rule_len = args.max_rule_len
alpha = args.alpha
graph_path = args.graph_path
thr = args.thr
top_k_rules = args.top_k_rules
y_value = [args.y_value]   # must be formatted as a list
output_file_name = args.output_file_name
calculus_mode = args.calculus_mode
parse_mode = args.parse_mode

data_colnames = [c for c in pd.read_csv(data_path).columns]
Y_variable = data_colnames[-1]
data_colnames = data_colnames[:-1]
X = data_colnames[split_column:]
if split_column == 0:
    Z = []
else:
    Z = data_colnames[:split_column]

df_graph = pd.read_table(graph_path)
G = nx.from_pandas_edgelist(df_graph, 'Gene1', 'Gene2', 'Score')

if calculus_mode:
    ALLSTAR_rule_calculus(data_path, X, Z, y_value, Y_variable, max_rule_len, alpha, G, top_k_rules, thr, output_file_name)
elif parse_mode:
    print(ALLSTAR_rule_parser(output_file_name))