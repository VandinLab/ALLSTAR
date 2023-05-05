import numpy as np
import pandas as pd
import argparse
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import graph_utils
from math import comb
from effect_calculator import EffectCalculator

# Analyzes the results of the synthetic experiments


external = 100
ds_sizes = [100,250,500,1000,5000,10000,25000]
ds_names = ["100","250","500","1000","5000","","25000"]

# synthetic exp 1
if True:
    print("*****")
    print("Experiment 1: random rules")
    print("*****")
    rel_path = "reports/g1"
    vs = ["naive","reliable","corrected"]

    cols = ["Wall_Time","Positive","Samples","Rule","Test_n","termination_condition"]
    tot_cnt_pos = {}
    for v in vs:
        tot_cnt_pos[v] = []
        print("-- VERSION =",v)
        df = pd.read_csv(rel_path+"/G1_"+v+"_res.csv", sep = "\t")
        for size in ds_sizes:
            cnt_tot = 0
            cnt_pos = 0
            for t_n in range(10):
                df_s = df[df["Samples"]==size]
                df_s = df_s[df_s["Test_n"]==t_n]
                cnt_tot += 1
                if np.sum(df_s["Positive"]) > 0:
                    cnt_pos += 1
            print("Ds size =",size,"---", cnt_pos,"/",cnt_tot)
            tot_cnt_pos[v].append(cnt_pos)
    print()


    ccs = ["tab:blue","tab:red","tab:green","tab:orange","tab:purple","tab:cyan"]   
    mks = ["o","v","<",">","^","*"]

    if True:
        fig, ax = plt.subplots()
        for v in vs:
            c = ccs[vs.index(v)]
            m = mks[vs.index(v)]
            ax.plot(ds_sizes,tot_cnt_pos[v],m+"-.",color = c, label="Version "+str(v))
            # ax.plot(ds_sizes,effs_g[v],m+":",color = c)
            # ax.set_yscale('log')

        # ax.hlines(0.8,xmin = min(ds_sizes),xmax = max(ds_sizes), label = "Real effect", color = "tab:gray")
        ax.set_xscale('log')
        # ax.set_ylim(-.05,1.05)
        ax.set_xticks(ds_sizes)
        ax.set_xticklabels(ds_names)
        ax.set_ylabel("Runs with False Positives")
        ax.set_xlabel("Dataset size")
        ax.legend()

        fig.savefig("images/G1_res.png")

matplotlib.rcParams.update({'font.size': 17})

external = 15
ds_sizes = [25,50,75,100,250,500,1000,5000,10000,25000]
max_lengths = [3,4,5,6,7]
mean_wt_g = {}
mean_eff_g = {}
mean_wt_no_g = {}
mean_eff_no_g = {}
#synthetic exp 2 with no plots
if True:
    print("*****")
    print("Experiment 2: one implanted rule. Test of correctness + shorter times + better bounds")
    print("*****")
    rel_path = "reports/g2"
    vs = ["","graph_"]

    cols = ["Wall_Time","Score","Rule","Samples","Max_Rule_Length","Test_n","termination_condition"]

    print("* Check number of FPs")
    for v in vs:
        print("-- VERSION =",v, "(empty = no graph)")
        df = pd.read_csv(rel_path+"/G2_"+v+"e"+str(external)+"_res.csv", sep = "\t")
        
        for size in ds_sizes:
            cnt_tot = 0
            cnt_pos = 0
            for t_n in range(10):
                df_s = df[df["Samples"]==size]
                df_s = df_s[df_s["Test_n"]==t_n]
                cnt_tot += 1

                fp = False
                for e,r in df_s[["Score","Rule"]].values:
                    if "Y \in [0]" in r:
                        fp = fp or e>0
                    elif "Y \in [1]" in r:
                        if e > 0:
                            fp = fp or not r in ["AKT1_somatic \in [1] & TP53_somatic \in [1] -> Y \in [1]","TP53_somatic \in [1] & AKT1_somatic \in [1] -> Y \in [1]"]

                cnt_pos += fp*1
            print("Ds size =",size,"---", cnt_pos,"/",cnt_tot)

    print("* Check number of FPs")
    for v in vs:
        wt  = {}
        eff = {}
        print("-- VERSION =",v, "(empty = no graph)")
        df = pd.read_csv(rel_path+"/G2_"+v+"e"+str(external)+"_res.csv", sep = "\t")
        for m_l in max_lengths:
            wt[m_l] = []
            eff[m_l]= []  
            print("-- MAX_LEN =",m_l)
            df_m = df[df["Max_Rule_Length"]==m_l]
            for size in ds_sizes:
                df_s = df_m[df_m["Samples"]==size]
                wt[m_l].append(max(0,np.mean(df_s["Wall_Time"])))
                df_s = df_s[df_s["Score"]>0]
                eff[m_l].append(max(0,np.mean(df_s["Score"])))
                print("Ds size =",size,
                    "--- Mean wall time {:1.2f}".format(np.mean(df_s["Wall_Time"])),
                    "--- Mean rule eff {:1.2f}".format(np.mean(df_s["Score"])),
                    "--- Rules found", len(df_s))
        if v == "":
            mean_wt_no_g = wt
            mean_eff_no_g = eff
        else:
            mean_wt_g = wt
            mean_eff_g = eff

    print()

#plot for walltime in exp 2
if True:
    # plots for exp 2
    ccs = ["tab:blue","tab:red","tab:green","tab:orange","tab:purple","tab:cyan"]
    mks = ["o","v","<",">","^","*"]
    ds_names = ["25","50","","100","250","500","1000","5000","","25000"]
    if True:
        fig, ax = plt.subplots(figsize=(10,6))
        for l in max_lengths:
            c = ccs[max_lengths.index(l)]
            m = mks[max_lengths.index(l)]
            ax.plot(ds_sizes,mean_wt_no_g[l],m+"-.",color = c, label="Max. len "+str(l))
            ax.plot(ds_sizes,mean_wt_g[l],m+":",color = c)
        ax.hlines(60*10,xmin = min(ds_sizes),xmax = max(ds_sizes), label = "Max wall time", color = "tab:gray")
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xticks(ds_sizes)
        ax.set_xticklabels(ds_names)
        ax.set_yticks([0.1,1,60,600])
        ax.set_yticklabels(["0.1 s","1 s","1 m","10 m"])
        ax.set_ylabel("Time (s)")
        ax.set_xlabel("Dataset size")
        ax.legend()
        fig.savefig("images/G2_Mean_wall_times.png")

# plot for rule effect in exp 2
if True:
    # plots for exp 2
    mean_eff_g = {}
    mean_eff_no_g = {}
    ccs = ["tab:blue","tab:red","tab:green","tab:orange","tab:purple","tab:cyan"]
    if True:
        data_path =  "inputs/synthetic/g2/g2_S25_t0_e15.csv"
        graph_path = 'PPI/FIsInGene_122921_with_annotations.txt'
        include_unmutated_treats = False

        data = pd.read_csv(data_path)
        data_colnames = data.columns.values.tolist()
        Y_variable = data_colnames[-1]
        # Drop last column as it is the outcome 
        data_colnames = data_colnames[:-1]
        # Isolate treatments
        Xs = data_colnames[7:]
        Zs = data_colnames[:7]

        # Calculate effects
        cfnds = []
        treats = []
        treats_vals = [[]]

        Y_classes = [[1]]
        Y_values = len(Y_classes)
        G = None
        if True:
            df_graph = pd.read_table(graph_path)
            G = nx.from_pandas_edgelist(df_graph, 'Gene1', 'Gene2', 'Score')
        ########### NB per i figli inserisco solo i valori in Xs.
        ########### se voglio quindi tenere o togliere qualcosa, conviene farlo PRIMA di dividere la matrice
        ########### nelle varie X Y e Z
        effs_no_g = {}
        effs_g = {}
        for mr in [3,4,5,6,7]:
            for gr in [False,True]:
                alpha_set= .05
                max_rule_len = mr
                use_graph_per_children = gr

                if True:
                    print("- Starting MHT correction calculus")
                    lxs = len(Xs)
                    if use_graph_per_children:
                        number_hyp_per_level = graph_utils.calculate_children_number_w_graph(Xs,G,Zs,data, max_rule_len-1)  # non cumulative number of distinct hypotheses per level if using the Protein-Protein graph
                    else:
                        number_hyp_per_level = [comb(lxs,i) for i in range(1,max_rule_len)]  # non cumulative number of distinct hypotheses per level if using the Protein-Protein graph
                    stop_level = max_rule_len-1
                    cumulative_number_of_Xs = sum(number_hyp_per_level[:stop_level])

                    tot_N_hyps = cumulative_number_of_Xs * Y_values * (len(Zs) + 1)  
                    
                    if include_unmutated_treats:
                        tot_N_hyps = tot_N_hyps * 2

                    alpha_set = alpha_set/tot_N_hyps
                    print("- Ended MHT correction calculus")

                means = []
                for d_s in ds_sizes:
                    effs = []
                    for t_n in range(10):
                        data_path =  "inputs/synthetic/g2/g2_S"+str(d_s)+"_t"+str(t_n)+"_e15.csv"
                        df = pd.read_csv(data_path)
                        for y_v in Y_classes:
                            eff = EffectCalculator.calculate_effect(df, [], ['TP53_somatic', 'AKT1_somatic', 'ANK2_somatic', 'TAF9_somatic', 'BRCA1_somatic', 'BRCA2_somatic'], [[1],[1],[1],[1],[1],[1]], Y_variable, [1], laplace_correction = True, return_per_strata_data = False, return_per_strata_a_counts = False,
                                            apply_reliable_correction = True, correction_alpha = alpha_set)
                            effs.append(eff)
                    means.append(np.mean(effs))

                if gr:
                    effs_g[mr] = means
                else:
                    effs_no_g[mr] = means
                    

        fig, ax = plt.subplots(figsize=(10,6))
        for l in max_lengths:
            c = ccs[max_lengths.index(l)]
            m = mks[max_lengths.index(l)]
            ax.plot(ds_sizes,effs_no_g[l],m+"-.",color = c, label="Max. len "+str(l))
            ax.plot(ds_sizes,effs_g[l],m+":",color = c)
            # ax.set_yscale('log')

        ax.hlines(0.88,xmin = min(ds_sizes),xmax = max(ds_sizes), label = "True effect", color = "tab:gray")
        ax.set_xscale('log')
        ax.set_ylim(-.05,1.05)
        ax.set_xticks(ds_sizes)
        ax.set_xticklabels(ds_names)
        ax.set_ylabel("Effect")
        ax.set_xlabel("Dataset size")
        ax.legend()

        fig.savefig("images/G2_Mean_effs_v2.png")


external = 100
ds_sizes = [100,250,500,1000,5000,10000,25000]
# synthetic exp 3 
if True:
    print("*****")
    print("Experiment 3: distinct rules")
    print("*****")
    rel_path = "reports/g3"
    vs = ["light","strong"]
    tot_cnt_pos = {}
    cols = ["Rule","Score","Parent","Level","Samples","Test_n","termination_condition"]
    print("3 distinct rules")
    for v in vs:
        tot_cnt_pos[v] = []
        print("-- VERSION =",v)
        df = pd.read_csv(rel_path+"/G3_"+v+"_res.csv", sep = "\t")
        for size in ds_sizes:
            cnt_tot = 0
            cnt_pos = 0
            for t_n in range(10):
                df_s = df[df["Samples"]==size]
                df_s = df_s[df_s["Test_n"]==t_n]
                cnt_tot += 1
                rules = [r.replace("_clone","") for r in df_s["Rule"].values]
                cnt_pos += 1*(len(rules)!=len(set(rules)))
            print("Ds size =",size,"---", cnt_pos,"/",cnt_tot)
            tot_cnt_pos[v].append(cnt_pos)
    print("Last rule useless cnt")
    for v in vs:
        print("-- VERSION =",v)
        df = pd.read_csv(rel_path+"/G3_"+v+"_res.csv", sep = "\t")
        for size in ds_sizes:
            cnt_tot = 0
            cnt_pos = 0
            for t_n in range(10):
                df_s = df[df["Samples"]==size]
                df_s = df_s[df_s["Test_n"]==t_n]
                df_s = df_s[df_s["Level"]==3]
                cnt_tot +=1
                cnt_pos += 1 * (df_s["Score"].values[0]<0)
            print("Ds size =",size,"---", cnt_pos,"/",cnt_tot)
    print()


    ds_names = ["100","250","500","1000","5000","10000","25000"]

    # plots synthetic exp 3
    if True:

        ccs = ["tab:blue","tab:red","tab:green","tab:orange","tab:purple","tab:cyan"]   
        mks = ["o","v","<",">","^","*"]

        if True:
            fig, ax = plt.subplots()
            for v in vs:
                c = ccs[vs.index(v)]
                m = mks[vs.index(v)]
                ax.plot(ds_sizes,tot_cnt_pos[v],m+"-.",color = c, label="Version "+str(v))
                # ax.plot(ds_sizes,effs_g[v],m+":",color = c)
                # ax.set_yscale('log')

            # ax.hlines(0.8,xmin = min(ds_sizes),xmax = max(ds_sizes), label = "Real effect", color = "tab:gray")
            ax.set_xscale('log')
            # ax.set_ylim(-.05,1.05)
            ax.set_xticks(ds_sizes)
            ax.set_xticklabels(ds_names)
            ax.set_ylabel("Runs with at least one clone returned")
            ax.set_xlabel("Dataset size")
            ax.legend()

            fig.savefig("images/G3_res.png")

