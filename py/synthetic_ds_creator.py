import argparse
import os
import pandas as pd
import numpy as np
from tqdm import tqdm

# Dataset creator for synthetic experiments

parser = argparse.ArgumentParser('main', description="Launch the causal rule discovery tool")
parser.add_argument('-v', '--version', type=int, required=True)
parser.add_argument('-r', '--repetitions', type=int, required=True)
parser.add_argument('-e', '--ext_vars', type=int, required=True)

args=parser.parse_args()

version = args.version
repetitions = args.repetitions
ext_vars = args.ext_vars

if not os.path.exists("inputs"):
    os.mkdir("inputs")
if not os.path.exists("inputs/synthetic/"):
    os.mkdir("inputs/synthetic/")
if not os.path.exists("inputs/synthetic/"+"g"+str(version)):
    os.mkdir("inputs/synthetic/"+"g"+str(version))
base_path = "inputs/synthetic/"+"g"+str(version)

if True:
    name_list = ['TP53_somatic',
    'AKT1_somatic',
    'ANK2_somatic',
    'TAF9_somatic',
    'BRCA1_somatic',
    'BRCA2_somatic',
    'CHD4_somatic',
    'CREBBP_somatic',
    'ERBB2_somatic',
    'FOXA1_somatic',
    'HTT_somatic',
    'MTOR_somatic',
    'MUC2_somatic',
    'PRKDC_somatic',
    'PTEN_somatic',
    'RB1_somatic',
    'RLF_somatic',
    'TAF1_somatic',
    'TAF1L_somatic',
    'TRRAP_somatic',
    'ZNF292_somatic',
    'TENM2_somatic',
    'HSPG2_somatic',
    'UNC13C_somatic',
    'SEMA3D_somatic',
    'ROS1_somatic',
    'LAMA2_somatic',
    'HERC2_somatic',
    'XIRP2_somatic',
    'TEP1_somatic',
    'SYNE1_somatic',
    'CNTNAP5_somatic',
    'LRP2_somatic',
    'ANKRD17_somatic',
    'CENPE_somatic',
    'MUC16_somatic',
    'TACC2_somatic',
    'KCNT2_somatic',
    'PLEC_somatic',
    'KIAA1210_somatic',
    'FASN_somatic',
    'PHKA2_somatic',
    'ANK1_somatic',
    'THSD7B_somatic',
    'ASPM_somatic',
    'MYLK_somatic',
    'PIK3CA_somatic',
    'DNAH2_somatic',
    'FREM2_somatic',
    'CSPP1_somatic',
    'DST_somatic',
    'PCLO_somatic',
    'PTPRD_somatic',
    'FAT2_somatic',
    'RNF213_somatic',
    'DMXL2_somatic',
    'VPS13D_somatic',
    'ASTN1_somatic',
    'SHANK1_somatic',
    'CMYA5_somatic',
    'COL14A1_somatic',
    'KMT2C_somatic',
    'NRXN2_somatic',
    'HUWE1_somatic',
    'ARID1B_somatic',
    'RELN_somatic',
    'RP1_somatic',
    'GON4L_somatic',
    'VWF_somatic',
    'CIT_somatic',
    'KDM6A_somatic',
    'SMG1_somatic',
    'MUC4_somatic',
    'TSHZ3_somatic',
    'ALMS1_somatic',
    'CDH1_somatic',
    'ABCA13_somatic',
    'MAP3K1_somatic',
    'ICE1_somatic',
    'AKAP9_somatic',
    'ZFHX4_somatic',
    'SSPOP_somatic',
    'FCGBP_somatic',
    'SVIL_somatic',
    'DNAH9_somatic',
    'TEX15_somatic',
    'HMCN1_somatic',
    'VPS13C_somatic',
    'ARHGAP31_somatic',
    'SPTB_somatic',
    'HRNR_somatic',
    'BSN_somatic',
    'ZNF462_somatic',
    'CPLANE1_somatic',
    'EPG5_somatic',
    'CUBN_somatic',
    'FMN2_somatic',
    'ZFPM2_somatic',
    'MYO18B_somatic',
    'MYO10_somatic',
    'TLN1_somatic',
    'FLNC_somatic',
    'NCOA6_somatic',
    'ANK3_somatic',
    'MICAL2_somatic',
    'SI_somatic',
    'DNAH17_somatic',
    'UBR5_somatic',
    'LAMA1_somatic',
    'MDN1_somatic']

ds_sizes = [100,250,500,1000,5000,10000,25000]
if version in [2,4]:
    ds_sizes=[25,50,75]+ds_sizes

# useful for testing no rules (= comparison test 1)
if version == 1:
    a = 0.4
    b = 0.7
    c = 0.2
    d = 0.5
    for n_test in tqdm(range(repetitions)):
        for tot_size in ds_sizes:
            file_path = os.path.join(base_path,"g"+str(version)+"_S"+str(tot_size)+"_t"+str(n_test)+"_e"+str(ext_vars)+".csv")
            if not os.path.exists(file_path):
                x1 = np.random.binomial(1, a, tot_size)
                x2 = np.random.binomial(1, b, tot_size)
                x_gen1 = np.bitwise_and(x1,x2)
                y = np.bitwise_or(x_gen1, np.random.binomial(1, c, tot_size))
                # data = [x1,x2]   # did not include them so there is no implanted rules by just watching the data
                data = []
                for i in range(ext_vars):
                    data += [(np.random.binomial(1, d, tot_size)).tolist()]
                data += [y]
                data = np.reshape(np.array(data).T,(-1,1+ext_vars))
                cols_names = name_list[:ext_vars] + ["Y"]
                df = pd.DataFrame(data, columns = cols_names)
                df.to_csv(file_path,index=False)

# useful for testing just one rule (= comparison test 2), where the rule has six element
# I'd also use it for testing overall performance of the algorithm in terms of
# shift from the real effect and influence of dataset size
if version == 2:
    a = 0.9
    d = 0.5
    for n_test in tqdm(range(repetitions)):
        for tot_size in ds_sizes:
            file_path = os.path.join(base_path,"g"+str(version)+"_S"+str(tot_size)+"_t"+str(n_test)+"_e"+str(ext_vars)+".csv")
            if not os.path.exists(file_path):
                x1 = np.random.binomial(1, a, tot_size)
                x2 = np.random.binomial(1, a, tot_size)
                x3 = np.random.binomial(1, a, tot_size)
                x4 = np.random.binomial(1, a, tot_size)
                x5 = np.random.binomial(1, a, tot_size)
                x6 = np.random.binomial(1, a, tot_size)
                y = np.bitwise_and(x1, x2) 
                y = np.bitwise_and(y, x3) 
                y = np.bitwise_and(y, x4) 
                y = np.bitwise_and(y, x5) 
                y = np.bitwise_and(y, x6) 
                y = np.bitwise_and(y, np.random.binomial(1, a, tot_size))
                data = [x1, x2, x3, x4, x5, x6]
                for i in range(ext_vars):
                    data += [(np.random.binomial(1, d, tot_size)).tolist()]
                data += [y]
                data = np.reshape(np.array(data).T,(-1,6+1+ext_vars))
                cols_names = name_list[:6]+name_list[-ext_vars:]+["Y"]
                df = pd.DataFrame(data, columns = cols_names)
                df.to_csv(file_path,index=False)

# useful for testing the removal of variables similar to the ones of the best rule
if version == 3:
    Xs = ["X1","X2","X3","X4","X5"]+["X1_c","X2_c","X3_c","X4_c","X5_c"]
    a = 0.5
    b = 0.4
    c = 0.7
    d = 0.65
    e = 0.15
    ch_1 = 0.025
    ch_2 = 0.975
    rnd = 0.5

    for n_test in tqdm(range(repetitions)):
        for tot_size in ds_sizes:
            file_path = os.path.join(base_path,"g"+str(version)+"_S"+str(tot_size)+"_t"+str(n_test)+"_e"+str(ext_vars)+".csv")
            if not os.path.exists(file_path) or True:
                x1 = np.random.binomial(1, a, tot_size)
                x2 = np.random.binomial(1, b, tot_size)
                x_gen1 = np.bitwise_and(x1,x2)
                x3 = np.random.binomial(1, c, tot_size)
                x4 = np.random.binomial(1, d, tot_size)
                x_gen2 = np.bitwise_and(x3,x4)
                x5 = np.random.binomial(1, e, tot_size)
                y = np.bitwise_or(x_gen1, x_gen2, x5)
                data = [x1,x2,x3,x4,x5]
                # this makes a vector similar to the original one if there is a majority of ones
                x1_c = np.bitwise_xor(x1, np.random.binomial(1, ch_1, tot_size))
                x2_c = np.bitwise_xor(x2, np.random.binomial(1, ch_1, tot_size))
                # this makes a vector similar to the opposite if there is a majority of ones
                x3_c = np.bitwise_xor(x3, np.random.binomial(1, ch_2, tot_size))
                x4_c = np.bitwise_xor(x4, np.random.binomial(1, ch_2, tot_size))
                x5_c = np.bitwise_and(x5, np.random.binomial(1, ch_2, tot_size))
                data += [x1_c,x2_c,x3_c,x4_c,x5_c]
                for i in range(ext_vars):
                    data += [(np.random.binomial(1, rnd, tot_size)).tolist()]
                data += [y]
                data = np.reshape(np.array(data).T,(-1,1+5*2+ext_vars))
                cols_names = name_list[:5]+[nn+"_clone" for nn in name_list[:5]]+name_list[5:ext_vars+5]+["Y"]
                df = pd.DataFrame(data, columns = cols_names)
                df.to_csv(file_path,index=False)

