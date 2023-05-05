import os
import pandas as pd
import sys
import re
import argparse
from effect_calculator import EffectCalculator

# Parser for synthetic experiments and save the report in reports/

parser=argparse.ArgumentParser('main', description="Parse the synthetic results")
parser.add_argument('-i', '--input_dir', required=True)
parser.add_argument('-v', '--version', type=int, required=True)
parser.add_argument('--naive', action='store_true')
parser.add_argument('--reliable', action='store_true')
parser.add_argument('--corrected', action='store_true')
parser.add_argument('--graph', action='store_true')
parser.add_argument('--no_graph', action='store_true')
parser.add_argument('-e','--external', type=int, default = -1)
parser.add_argument('--light', action='store_true')
parser.add_argument('--strong', action='store_true')

args=parser.parse_args()

version = args.version
input_dir = args.input_dir

naive = args.naive
reliable = args.reliable
corrected = args.corrected

graph = args.graph
no_graph = args.no_graph
external = args.external

light = args.light
strong = args.strong


if not os.path.exists("reports"):
    os.mkdir("reports")
if not os.path.exists("reports/g"+str(version)):
    os.mkdir("reports/g"+str(version))
output_dir = "reports/g"+str(version)

assert(version not in [2] or external!=-1)

if version == 1:
    results = pd.DataFrame(columns = ['Wall_Time', 'Positive','Samples'])

    files = [f for f in os.listdir(input_dir)]

    # select only files of the runs of interest
    if naive:
        files = [f for f in files if "_naive" in f]
    if reliable:
        files = [f for f in files if "_reliable" in f]
    if corrected:
        files = [f for f in files if "_corrected" in f]


    for filename in files:
        n_samples = re.search('g1_S(.*)_t', filename)
        positive = 0
        with open(os.path.join(input_dir, filename), 'r') as f:
            data = f.readlines()
            for line in data:
                if '- objective:' in line:
                    obj = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)
                    if float(obj[0]) > 0:
                        positive = 1
                        #print(float(obj[0]))
                if '- wall_time' in line:
                    wall_time = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)
                if '- termination_condition' in line:
                    termination_condition = line.split(":")[1]
                if "('" in line[:5]:                
                    res = eval(line)
                    rule = EffectCalculator.rule_tostring(cfnds=res[1], treats=res[2], treats_vals=res[3],
                            outcome=res[4],outcome_vals=res[5])

                    results = results.append({
                                        'Rule' : rule, 
                                        'Wall_Time' : float(wall_time[0]), 
                                        'Positive' : str(positive),
                                        'Samples': str(n_samples.group(1)),
                                        'Test_n': filename.split("_t")[1].split("_")[0],
                                        'termination_condition':termination_condition
                                        }, ignore_index = True)

    res_name = "G1"
    if naive:
        res_name += "_naive"
    if reliable:
        res_name += "_reliable"
    if corrected:
        res_name += "_corrected"

    results.to_csv(output_dir + '/' + res_name + '_res.csv', sep ='\t')

if version == 2:
    results = pd.DataFrame(columns = ['Wall_Time', 'Score','Rule','Samples','Test_n','Max_Rule_Length'])

    files = [f for f in os.listdir(input_dir) if "_e"+str(external) in f]

    # select only files of the runs of interest
    if graph:
        files = [f for f in files if "_graph" in f]
    if no_graph:
        files = [f for f in files if "_graph" not in f]


    for filename in files:
        n_samples = re.search('g2_S(.*)_t', filename)
        if 'max_len_3' in filename:
            mrl = 3
        if 'max_len_4' in filename:
            mrl = 4
        if 'max_len_5' in filename:
            mrl = 5
        if 'max_len_6' in filename:
            mrl = 6
        if 'max_len_7' in filename:
            mrl = 7
        with open(os.path.join(input_dir, filename), 'r') as f:
            data = f.readlines()
            for line in data:
                if '- objective:' in line:
                    obj = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)
                if '- wall_time' in line:
                    wall_time = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)
                    if "ms" in line:
                        div = 1000
                    elif " m" in line:
                        div = 1/60
                    elif " s" in line:
                        div = 1
                if '- termination_condition' in line:
                    termination_condition = line.split(":")[1]

                if "('" in line[:5]:
                    res = eval(line)
                    rule = EffectCalculator.rule_tostring(cfnds=res[1], treats=res[2], treats_vals=res[3],
                            outcome=res[4],outcome_vals=res[5])
                    results = results.append({
                                        'Wall_Time' : float(wall_time[0])/div, 
                                        'Score' : str(float(obj[0])),
                                        'Rule': rule,
                                        'Samples': str(n_samples.group(1)),
                                        'Test_n': filename.split("_t")[1].split("_")[0],
                                        'Max_Rule_Length': str(mrl),
                                        'termination_condition':termination_condition
                                        }, ignore_index = True)


    res_name = "G2"
    if graph:
        res_name += "_graph"

    results.to_csv(output_dir + '/' + res_name +"_e"+str(external) + '_res.csv', sep ='\t')

if version == 3:
    results = pd.DataFrame(columns = ['Rule', 'Score','Parent','Level', 'Samples'])
    directory_name_rules = "results/synth_v13/v13/"
    files = [f for f in os.listdir(directory_name_rules)]

    if light:
        files = [f for f in files if "_light" in f]
    if strong:
        files = [f for f in files if "_strong" in f]

    for filename in files:
        parent_rule = ""
        eff = 0
        lv = 0
        with open(os.path.join(directory_name_rules, filename), 'r') as f:
            data = f.readlines()
            samples = int(filename.split("_S")[1].split("_")[0])
            for line in data:
                if '- objective:' in line:
                    obj = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)
                if '- termination_condition' in line:
                    termination_condition = line.split(":")[1]

                if "('" in line[:5]:
                    res = eval(line)
                    rule = EffectCalculator.rule_tostring(cfnds=res[1], treats=res[2], treats_vals=res[3],
                            outcome=res[4],outcome_vals=res[5])
                    if parent_rule== "":
                        parent_rule = rule
                    eff = float(obj[0])
                    results = results.append({
                                        'Rule' : rule, 
                                        'Score' : float(obj[0]),
                                        'Parent': parent_rule,
                                        'Level': lv,
                                        'Test_n': filename.split("_t")[1].split("_")[0],
                                        'Samples' : samples,
                                        'termination_condition':termination_condition
                                        }, ignore_index = True)
                    lv += 1

    res_name = "G3"
    if light:
        res_name += "_light"
    if strong:
        res_name += "_strong"

    results.to_csv(output_dir + '/' + res_name + '_res.csv', sep ='\t')
