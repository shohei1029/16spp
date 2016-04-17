#1/usr/bin/env python3

# 2016.4.17
# 配列領域ごとのDmeanを計算してplotする。
# Dmeanの計算は，別スクリプトであるcal_Dmean_by_sim.pyをimportして使っている。

import argparse

import pandas as pd
import matplotlib #to set use('Agg') 
matplotlib.use('Agg') #place this before any other pylab/matplotlib/pyplot import
import matplotlib.pyplot as plt
import seaborn as sns

import cal_Dmean_by_sim

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--target_seqregion", type=str, help="target sequence region. sep=,")
parser.add_argument("-o", "--output_file", type=str, default="./out_plot_bar.png"))
args = parser.parse_args()


def plot_barplot(df, outfile='./out_plot_bar.png'):
#    sns.set(style="whitegrid")
    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    fig, ax = plt.subplots(1, 1, figsize=(8,6)) #fig->figure obj. ax->graph obj. 2,1とかだとgは配列に.2,2だとarray.
    ax = sns.barplot(x='seq_region', y='Dmean', data=df)
    fig.savefig(outfile, dpi=180)
    sns.plt.close()


if __name__ == '__main__':
    if args.target_seqregion:
        seqregions = args.target_seqregion.split(',')
    else:
        seqregions = "RVP RVT_1 RVT_thumb RVT_connect RNase_H Integrase_Zn rve IN_DBD_C".split(' ')

    dmeans = []

    in_dir = "../analysis/sim/pol-4secs/"
    file_name_core = "_sec_A,B,C_blastp_1e-5_pol-4secs_A,B,C_HIV-1-gM-noRs_pol-aa_v3.sim.txt"

    for seqre in seqregions:
        in_file = in_dir + seqre + file_name_core

        print(in_file)
        in_fh = open(in_file, 'r')
        dmean = cal_Dmean_by_sim.cal_Dmean_file(in_fh)
        print(dmean)
        dmeans.append(dmean)

    df = pd.DataFrame({'seq_region':seqregions, 'Dmean':dmeans})
    plot_barplot(df, args.output_file)
