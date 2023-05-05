import numpy as np
import pandas as pd
from typing import Dict, List
import matplotlib.pyplot as plt
import seaborn as sb
import re
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")



results = pd.DataFrame(columns = ['CT_pval','rank_pval','CT_odds','rank_odds','CT_eff','rank_eff'])

ranking_file = 'results/brute_force_ranks_fisher_exact.log'
with open(ranking_file, 'r') as f:
        data = f.readlines()[6:]
        for line in tqdm(data):
            obj = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)
            results = results.append({
                                    'CT_pval' : str('['+obj[2]+','+obj[3]+','+obj[4]+','+obj[5]+']'),
                                    'rank_pval' : int(obj[1]),
                                    'CT_odds' : str('['+obj[8]+','+obj[9]+','+obj[10]+','+obj[11]+']'),
                                    'rank_odds' : int(obj[7]),
                                    'CT_eff' : str('['+obj[14]+','+obj[15]+','+obj[16]+','+obj[17]+']'),
                                    'rank_eff' : int(obj[13])
                                    }, ignore_index = True)

results.to_csv('results/parsed_rankings.csv', sep ='\t')


rankings = pd.read_csv('results/parsed_rankings.csv', sep="\t")
rankings.head()

r_pval = rankings[['CT_pval','rank_pval']]
r_odds = rankings[['CT_odds', 'rank_odds']]
r_eff = rankings[['CT_eff', 'rank_eff']]

rankings_1 = pd.merge(r_pval, r_odds, left_on='CT_pval', right_on='CT_odds', how='inner')
rankings_merged = pd.merge(rankings_1, r_eff, left_on='CT_pval', right_on='CT_eff', how='inner')
rankings_merged.drop(['CT_odds', 'CT_eff'], axis=1, inplace=True)
rankings_merged.rename(columns={'CT_pval': 'contingency'}, inplace=True)
rankings_merged.sort_values(by='rank_eff').head(20)

head_filter = 10
cmap = sb.diverging_palette(0, 230, 90, 60, as_cmap=True)

# Rank by effect
rankings_merged.sort_values(by='rank_eff', inplace=True)
fig, ax = plt.subplots(figsize=(11, 9))
sb.heatmap(rankings_merged.head(head_filter).drop('contingency', axis=1), cmap=cmap)
# xticks
yticks_labels = rankings_merged['contingency'].head(head_filter)
plt.yticks(np.arange(len(yticks_labels)) + .5, labels=yticks_labels, rotation='horizontal')
#plt.title('Sorted by effect: rank comparison', loc='center')
fig.savefig("images/fisher_ratios_effects_comparison.png")