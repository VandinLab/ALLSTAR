import statistics
import numpy as np
from collections import Counter
import pandas as pd
import math
import scipy.stats as st
import warnings

warnings.filterwarnings("ignore")

# Library used to calculate rule effects and contingency tables

class EffectCalculator:

    @staticmethod
    def calculate_effect(df, cfnds, treats, treats_vals, outcome, outcome_vals, laplace_correction = True, return_per_strata_data = False, return_per_strata_a_counts = False,
                            apply_reliable_correction = False, correction_alpha = 0.05):

        if len(treats) != len(treats_vals):
            print(treats,treats_vals)
        assert len(treats) == len(treats_vals)

        effect = 0
        tot_cnt = 0
        if cfnds == []:
            cnt_cfnds = {(0,):len(df)}
        else:
            cnt_cfnds = Counter(list(df[cfnds].itertuples(index=False, name=None)))

        per_strata_a_values = []
        per_strata_data = []

        for cf_vals in cnt_cfnds.keys():

            stratum_df = df.copy()
            # filter out the df to select only rows with specific Z values
            for i in range(len(cfnds)):
                stratum_df = stratum_df[stratum_df[cfnds[i]] == cf_vals[i]]

            stratum_df = stratum_df[treats+[outcome]]

            pos_rule = stratum_df.copy()
            for i in range(len(treats)):
                pos_rule = pos_rule[pos_rule[treats[i]].isin(treats_vals[i])]
            cnt_out = Counter(pos_rule[outcome].isin(outcome_vals))
            a = cnt_out[True]
            b = cnt_out[False]
            neg_rule = stratum_df.drop(stratum_df[stratum_df.index.isin(pos_rule.index)].index)
            cnt_out = Counter(neg_rule[outcome].isin(outcome_vals))
            c = cnt_out[True]
            d = cnt_out[False]

            # NB: data is the unbiased one, that is without correction
            if return_per_strata_a_counts:
                per_strata_a_values.append(a)
            if return_per_strata_data:
                per_strata_data.append((a,b,c,d,cnt_cfnds[cf_vals]))
            if laplace_correction:
                a += 1
                b += 1
                c += 1
                d += 1
                cnt_cfnds[cf_vals] += 4
            tot_cnt += cnt_cfnds[cf_vals]
            # for avoiding division by 0 errors if laplace = False
            if a+b == 0:
                b+=1
            if c+d == 0:
                d+=1
            strata_effect = a/(a+b) - c/(c+d) 
            if apply_reliable_correction:
                strata_effect -= EffectCalculator.apply_correction(correction_alpha, 2*math.sqrt(a+b), 2*math.sqrt(c+d))
            effect += strata_effect*cnt_cfnds[cf_vals]

        effect = effect/tot_cnt
        if return_per_strata_a_counts:
            return effect,per_strata_a_values
        if return_per_strata_data:
            return effect,per_strata_data
        return effect

    @staticmethod
    def apply_correction(correction_alpha, den1, den2):
        beta = st.norm.ppf(1-correction_alpha/2)
        return beta/den1 + beta/den2

    @staticmethod
    def rule_tostring(cfnds, treats, treats_vals, outcome, outcome_vals):
        s = ""
        for i in range(len(treats)):
            if s!="":
                s+=" & "
            s+= treats[i]+" \in "+str(treats_vals[i])
        s+= " -> "+outcome+" \in "+str(outcome_vals)
        if cfnds != []:
            s += " | "+"".join([str(e)+"," for e in cfnds])
            s = s[:-1]
        return s

    @staticmethod
    def get_z_strata_tot_values(df, cfnds, treats, treats_vals, outcome, outcome_vals):

        if len(treats) != len(treats_vals):
            print(treats,treats_vals)
        assert len(treats) == len(treats_vals)

        if cfnds == []:
            cnt_cfnds = {(0,):len(df)}
        else:
            cnt_cfnds = Counter(list(df[cfnds].itertuples(index=False, name=None)))
        per_strata_a_values = []

        for cf_vals in cnt_cfnds.keys():
            stratum_df = df.copy()
            # filter out the df to select only rows with specific Z values
            for i in range(len(cfnds)):
                stratum_df = stratum_df[stratum_df[cfnds[i]] == cf_vals[i]]

            stratum_df = stratum_df[treats+[outcome]]

            pos_rule = stratum_df.copy()
            for i in range(len(treats)):
                pos_rule = pos_rule[pos_rule[treats[i]].isin(treats_vals[i])]
            cnt_out = Counter(pos_rule[outcome].isin(outcome_vals))
            a = cnt_out[True]
            b = cnt_out[False]
            neg_rule = stratum_df.drop(stratum_df[stratum_df.index.isin(pos_rule.index)].index)
            cnt_out = Counter(neg_rule[outcome].isin(outcome_vals))
            c = cnt_out[True]
            d = cnt_out[False]
            # NB: data is the unbiased one, that is without correction
            per_strata_a_values.append((a+c,a+b+c+d))
            
        return per_strata_a_values

def tester():
    df = pd.read_csv("csvs/paper_example_rebranded.csv")

    cfnds = []
    treats = ["X"]
    treats_vals = [[2,3]]
    outcome = "Y"
    outcome_vals = [5]
    ss = EffectCalculator.rule_tostring(cfnds, treats, treats_vals, outcome, outcome_vals)
    v = EffectCalculator.calculate_effect(df, cfnds, treats, treats_vals, outcome, outcome_vals, laplace_correction = False)
    print(ss)
    print(v,(df[outcome]-4).mean())

    cfnds = []
    treats = ["E1"]
    treats_vals = [[0,1]]
    outcome = "Y"
    outcome_vals = [5]
    ss = EffectCalculator.rule_tostring(cfnds, treats, treats_vals, outcome, outcome_vals)
    v = EffectCalculator.calculate_effect(df, cfnds, treats, treats_vals, outcome, outcome_vals, laplace_correction = False)
    print(ss)
    print(v,(df[outcome]-4).mean())

    cfnds = ["Z"]
    treats = ["X"]
    treats_vals = [[2,3]]
    outcome = "Y"
    outcome_vals = [5]
    ss = EffectCalculator.rule_tostring(cfnds, treats, treats_vals, outcome, outcome_vals)
    v = EffectCalculator.calculate_effect(df, cfnds, treats, treats_vals, outcome, outcome_vals, laplace_correction = False)
    print(ss)
    print(v,(df[outcome]-4).mean())

    cfnds = ["Z"]
    treats = ["X"]
    treats_vals = [[3]]
    outcome = "Y"
    outcome_vals = [5]
    ss = EffectCalculator.rule_tostring(cfnds, treats, treats_vals, outcome, outcome_vals)
    v = EffectCalculator.calculate_effect(df, cfnds, treats, treats_vals, outcome, outcome_vals, laplace_correction = False)
    ss2 = EffectCalculator.rule_tostring(cfnds+["E1"], treats, treats_vals, outcome, outcome_vals)
    v2 = EffectCalculator.calculate_effect(df, cfnds+["E1"], treats, treats_vals, outcome, outcome_vals, laplace_correction = False)
    print(ss,"    vs    ",ss2)
    print(v,v2,0.18) # expected value for the rule is 0.18

    cfnds = ["Z"]
    treats = ["X"]
    treats_vals = [[3]]
    outcome = "Y"
    outcome_vals = [5]
    ss = EffectCalculator.rule_tostring(cfnds, treats, treats_vals, outcome, outcome_vals)
    v = EffectCalculator.calculate_effect(df, cfnds, treats, treats_vals, outcome, outcome_vals, laplace_correction = False)
    ss2 = EffectCalculator.rule_tostring(cfnds, treats, [[2]], outcome, outcome_vals)
    v2 = EffectCalculator.calculate_effect(df, cfnds, treats, [[2]], outcome, outcome_vals, laplace_correction = False)
    print(ss,"    vs    ",ss2)
    print(v,v2) 

    cfnds = ["Z"]
    treats = ["X"]
    treats_vals = [[3]]
    outcome = "Y"
    outcome_vals = [4,5]
    ss = EffectCalculator.rule_tostring(cfnds, treats, treats_vals, outcome, outcome_vals)
    v = EffectCalculator.calculate_effect(df, cfnds, treats, treats_vals, outcome, outcome_vals, laplace_correction = False)
    print(ss)
    print(v) 

    cfnds = ["Z"]
    treats = ["X"]
    treats_vals = [[3]]
    outcome = "Y"
    outcome_vals = [-4,-5]
    ss = EffectCalculator.rule_tostring(cfnds, treats, treats_vals, outcome, outcome_vals)
    v = EffectCalculator.calculate_effect(df, cfnds, treats, treats_vals, outcome, outcome_vals, laplace_correction = False)
    print(ss)
    print(v) 

    cfnds = ["Z"]
    treats = ["X"]
    treats_vals = [[3]]
    outcome = "Y"
    outcome_vals = [5]
    ss = EffectCalculator.rule_tostring(cfnds, treats, treats_vals, outcome, outcome_vals)
    v = EffectCalculator.calculate_effect(df, cfnds, treats, treats_vals, outcome, outcome_vals, laplace_correction = False)
    ss2 = EffectCalculator.rule_tostring(cfnds, treats+["E1"], treats_vals+[[0]], outcome, outcome_vals)
    v2 = EffectCalculator.calculate_effect(df, cfnds, treats+["E1"], treats_vals+[[0]], outcome, outcome_vals, laplace_correction = False)
    print(ss,"    vs    ",ss2)
    print(v,v2) 

    cfnds = ["Z"]
    treats = ["E1"]
    treats_vals = [[0]]
    outcome = "Y"
    outcome_vals = [5]
    ss = EffectCalculator.rule_tostring(cfnds, treats, treats_vals, outcome, outcome_vals)
    v = EffectCalculator.calculate_effect(df, cfnds, treats, treats_vals, outcome, outcome_vals, laplace_correction = False)
    print(ss)
    print(v) 


if __name__ == "__main__":
    tester()