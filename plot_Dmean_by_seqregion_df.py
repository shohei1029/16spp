#1/usr/bin/env python3

# 2016.4.17
# ./plot_Dmean_by_seqregion.py
# 配列領域ごとのDmeanを計算してplotする。
# Dmeanの計算は，別スクリプトであるcal_Dmean_by_sim.pyをimportして使っている。

# 2016.4.19
# plot_Dmean_by_seqregion_df.py へcp。cal_Dmean_by_sim.pyのimportをやめる。
# 配列間距離を全部dfの表としてまとめる。(消費メモリやばそう。。)。
# Dmean 計算もseaborn側にさせる
# ->すべてのエッジがsimファイルに書かれている前提。エッジカットがあると平均がちゃんと全配列間の平均にならない。

# 2016.4.25
# hard coding部分をなくし(argsのdefaultへ)，env等他タンパク質でも使えるように。


import os
import argparse

import pandas as pd
import matplotlib #to set use('Agg') 
matplotlib.use('Agg') #place this before any other pylab/matplotlib/pyplot import
import matplotlib.pyplot as plt
import seaborn as sns

import cal_Dmean_by_sim

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--target_seqregion", type=str, help="target sequence region. sep=,", default="RVP RVT_1 RVT_thumb RVT_connect RNase_H Integrase_Zn rve IN_DBD_C")
parser.add_argument("-o", "--output_dir", type=str, default="../results/Dmean/pol-domains-HXB2-Pfam")
parser.add_argument("-i", "--input_dir", type=str, default="../analysis/sim/pol-domains-HXB2-Pfam")
parser.add_argument("-n", "--name_core", type=str, default="_pi0_blastp_1e-5_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.txt")
args = parser.parse_args()


def parse_edge_weight(in_fh):
    dists = []
    ids = set()
    for line in in_fh:
        line = line.rstrip()
        id1, id2, sim = line.split(' ')
        dissim = 1 - float(sim)
        dists.append(dissim)
    
    return dists

def plot_Dmean_violinplot(df, outfile='./out_violinplot_Dmean.png'):
#    sns.set(style="whitegrid")
    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    fig, ax = plt.subplots(1, 1, figsize=(8,6)) #fig->figure obj. ax->graph obj. 2,1とかだとgは配列に.2,2だとarray.
    ax = sns.violinplot(x='seq_region', y='distance', data=df)
    plt.xlabel('sequence region')
    plt.ylabel('Dmean')
    fig.savefig(outfile, dpi=180)
    sns.plt.close()

def plot_Dmean_boxplot(df, outfile='./out_boxplot_Dmean.png'):
#    sns.set(style="whitegrid")
    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    fig, ax = plt.subplots(1, 1, figsize=(8,6)) #fig->figure obj. ax->graph obj. 2,1とかだとgは配列に.2,2だとarray.
    ax = sns.boxplot(x='seq_region', y='distance', data=df)
    plt.xlabel('sequence region')
    plt.ylabel('Dmean')
    fig.savefig(outfile, dpi=180)
    sns.plt.close()

def plot_Dmean_barplot(df, outfile='./out_barplot_Dmean.png'):
#    sns.set(style="whitegrid")
    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    fig, ax = plt.subplots(1, 1, figsize=(8,6)) #fig->figure obj. ax->graph obj. 2,1とかだとgは配列に.2,2だとarray.
    ax = sns.barplot(x='seq_region', y='distance', data=df)
    plt.xlabel('sequence region')
    plt.ylabel('Dmean')
    fig.savefig(outfile, dpi=180)
    sns.plt.close()


if __name__ == '__main__':
    seqregions = args.target_seqregion.split(',')

    outdir = args.output_dir
    os.makedirs(outdir, exist_ok=True)

    out_file_violin = outdir + "/out_violinplot_Dmean.png"
    out_file_box    = outdir + "/out_boxplot_Dmean.png"
    out_file_bar    = outdir + "/out_barplot_Dmean.png"

    in_dir = args.input_dir
    file_name_core = args.name_core

    dmeans = []
    variances = []
    for seqre in seqregions:
        in_file = in_dir + seqre + file_name_core

        print(in_file)
        in_fh = open(in_file, 'r')
        dists = parse_edge_weight(in_fh)

        try:
            df = df.append(pd.DataFrame({'seq_region':[seqre]*len(dists), 'distance':dists}))
        except NameError:
            df = pd.DataFrame({'seq_region':[seqre]*len(dists), 'distance':dists})

    plot_Dmean_violinplot(df, out_file_violin)
    plot_Dmean_boxplot(df, out_file_box)
    plot_Dmean_barplot(df, out_file_bar)


