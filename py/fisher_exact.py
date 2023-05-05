import math
import numpy as np
from typing import Dict, List
import argparse
import matplotlib.pyplot as plt
import scipy.stats as ss
from tqdm import tqdm
import os

# script to show that our effect calculus is not equivalent to standard Fisher or odds-based approaches

parser=argparse.ArgumentParser('main', description="Launch the causal rule discovery tool")
parser.add_argument('--graph_path', default = 'PPI/FIsInGene_122921_with_annotations.txt')
parser.add_argument('--use_MHT_correction', action='store_true')
parser.add_argument('--verbose', action='store_true')

args=parser.parse_args()

use_MHT_correction = args.use_MHT_correction
graph_path = args.graph_path
verbose = args.verbose

def plot_rankings(r1,r2,name_fig, xlabel, ylabel):
    
    fig, ax = plt.subplots()
    ax.plot(r1,r2,linestyle = 'None', marker = "o",markersize=1)

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    fig.savefig(name_fig)

if not os.path.exists("images"):
    os.mkdir("images")
if not os.path.exists("images/CMH_tests"):
    os.mkdir("images/CMH_tests")
if not os.path.exists("images/CMH_tests/ratios"):
    os.mkdir("images/CMH_tests/ratios")


# brute force true force
def calculate_abs_rule_val(a,b,c,d):
    beta = ss.norm.ppf(1-0.05/2)
    return np.abs((a+1)/(a+b+2)-(c+1)/(c+d+2))-beta/2/math.sqrt(a+b+2)-beta/2/math.sqrt(c+d+2)

pvs      = []
ratios_c = []
effs     = []

max_size = 100
for n_sigma in tqdm(range(max_size)):
    n_n_sigma = max_size - n_sigma
    for a in range(n_sigma):
        b = n_sigma - a
        for c in range(n_n_sigma):
            d = n_n_sigma - c

            ar = [[a+1,b+1],[c+1,d+1]]
            p_vs_FIS_scipy = min(ss.fisher_exact(ar,alternative="less")[1],
                                    ss.fisher_exact(ar,alternative="greater")[1])
            pvs.append((str(ar),p_vs_FIS_scipy))
            effs.append((str(ar),calculate_abs_rule_val(a,b,c,d)))
            ratios_c.append((str(ar),max((a+1)*(d+1)/(b+1)/(c+1),1/((a+1)*(d+1)/(b+1)/(c+1)))))
print("Generated")
only_pvs = [e[1] for e in pvs]
only_ratios_c = [e[1] for e in ratios_c]
only_effs = [e[1] for e in effs]

pvs_sorted = sorted(only_pvs)
ratios_c_sorted = sorted(only_ratios_c, reverse = True)
effs_sorted = sorted(only_effs, reverse = True)

ranks_pvs = [pvs_sorted.index(pvs[i][1]) for i in range(len(pvs))]
ranks_ratios_c = [ratios_c_sorted.index(ratios_c[i][1]) for i in range(len(ratios_c))]
ranks_effs = [effs_sorted.index(effs[i][1]) for i in range(len(effs))]

print("Ranked")

pvs = sorted(pvs,key=lambda x:(x[1],x[0]))
ratios_c = sorted(ratios_c,key=lambda x:(-x[1],x[0]), reverse = False) # in this way I preserve lexicographic order
effs = sorted(effs,key=lambda x:(-x[1],x[0]), reverse = False) # in this way I preserve lexicographic order
print("Sorted")
tables_pvs = [e[0] for e in pvs]
tables_ratios_c = [e[0] for e in ratios_c]
tables_effs = [e[0] for e in effs]

print("Kendall's tau for Fisher p values and effects",ss.kendalltau(ranks_pvs, ranks_effs))
plot_rankings(ranks_pvs, ranks_effs,"images/CMH_tests/ratios/ranks_FIS_effs.png", "Rank of Fisher p-values", "Rank of rule effects")
print("Kendall's tau for Fisher p values and ratios ",ss.kendalltau(ranks_pvs, ranks_ratios_c))
plot_rankings(ranks_pvs, ranks_ratios_c,"images/CMH_tests/ratios/ranks_FIS_odd_ratios.png", "Rank of Fisher p-values", "Rank of odds ratios")
print("Kendall's tau for odds ratios and effects    ",ss.kendalltau(ranks_ratios_c, ranks_effs))
plot_rankings(ranks_ratios_c, ranks_effs,"images/CMH_tests/ratios/ranks_odd_ratios_effs.png", "Rank of odds ratios", "Rank of rule effects")

for i in range(len(tables_pvs)):
    if len(set([tables_pvs[i],tables_ratios_c[i],tables_effs[i]]))!=1:
        print("line "+str(i+1)+":",pvs_sorted.index(pvs[i][1]),"-",tables_pvs[i],": pv=",pvs[i][1],"  ",ratios_c_sorted.index(ratios_c[i][1]),"-",tables_ratios_c[i],": ratio=",ratios_c[i][1],"  ",effs_sorted.index(effs[i][1]),"-",tables_effs[i],": eff=",effs[i][1])

