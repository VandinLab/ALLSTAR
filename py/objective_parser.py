import os
import pandas as pd
import sys
import re
import argparse
from effect_calculator import EffectCalculator

# Script to parse the results and create the aggregated report of multiple runs

parser=argparse.ArgumentParser('main', description="Parse the results")
parser.add_argument('-f', '--data_path', required=True)
parser.add_argument('-v', '--version', type=int, required=True)
parser.add_argument('--max_rule_len', type=int, default=4)
parser.add_argument('--use_MHT_correction', action='store_true')
parser.add_argument('--use_graph_per_children', action='store_true')
parser.add_argument('--include_unmutated_treats', action='store_true')

args=parser.parse_args()

data_path = args.data_path
max_rule_len = args.max_rule_len
use_MHT_correction = args.use_MHT_correction
include_unmutated_treats = args.include_unmutated_treats
use_graph_per_children = args.use_graph_per_children
version = args.version


if not os.path.exists("reports"):
    os.mkdir("reports")
dataset_name_only = data_path.split("/")[-1].split(".")[0]
if not os.path.exists("reports/"+dataset_name_only):
    os.mkdir("reports/"+dataset_name_only)
if not os.path.exists("reports/"+dataset_name_only+"/v"+str(version)):
    os.mkdir("reports/"+dataset_name_only+"/v"+str(version))
output_dir = "reports/"+dataset_name_only+"/v"+str(version)

results = pd.DataFrame(columns = ['Rule', 'Score','Parent','Level','Per Strata Table (a,b,c,d,n)'])
df_in = pd.read_csv(data_path)
dataset_name_only = data_path.split("/")[-1].split(".")[0]
directory_name_report =  "reports/"+dataset_name_only+"/v"+str(version)
directory_name_rules =  "results/"+dataset_name_only+"/v"+str(version)

files = [f for f in os.listdir(directory_name_rules) if "_max_len_"+str(max_rule_len) in f]

# select only files of the runs of interest
if include_unmutated_treats:
    files = [f for f in files if "_w_0s" in f]
else:
    files = [f for f in files if "_w_0s" not in f]
if use_MHT_correction:
    files = [f for f in files if "_corr" in f]
else:
    files = [f for f in files if "_corr" not in f]
if use_graph_per_children:
    files = [f for f in files if "_graph" in f]
else:
    files = [f for f in files if "_graph" not in f]

for filename in files:
    parent_rule = ""
    eff = 0
    lv = 0
    with open(os.path.join(directory_name_rules, filename), 'r') as f:
        data = f.readlines()
        for line in data:
            if eff < -0.2:
                break
            if '- objective:' in line:
                obj = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)
            if "('" in line[:5]:
                res = eval(line)
                rule = EffectCalculator.rule_tostring(cfnds=res[1], treats=res[2], treats_vals=res[3],
                        outcome=res[4],outcome_vals=res[5])
                if parent_rule== "":
                    parent_rule = rule
                effect,per_strata_data = EffectCalculator.calculate_effect(df_in, cfnds=res[1], treats=res[2], treats_vals=res[3], \
                    outcome=res[4],outcome_vals=res[5], laplace_correction = True, return_per_strata_data = True) # we don't care abt other effect parameters since we care only abt stratas that are fixed in in_df independently of corrections or not
                eff = float(obj[0])
                results = results.append({
                                    'Rule' : rule, 
                                    'Score' : float(obj[0]),
                                    'Parent': parent_rule,
                                    'Level': lv,
                                    'Per Strata Table (a,b,c,d,n)' : per_strata_data
                                    }, ignore_index = True)
                lv += 1

if max(results['Level']) == 0:  # meaning I have no "consecutive" cleaned rules
    results = results.sort_values("Score", ascending=False)

res_name = dataset_name_only+"_v"+str(version)
res_name += "_max_len_"+str(max_rule_len)
if include_unmutated_treats:
    res_name += "_w_0s"
if use_MHT_correction:
    res_name += "_corr"
if use_graph_per_children:
    res_name += "_graph"

results.to_csv(directory_name_report + '/' + res_name + '_res.csv', sep ='\t')