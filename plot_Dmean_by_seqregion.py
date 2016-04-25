#1/usr/bin/env python3

# 2016.4.17
# 配列領域ごとのDmeanを計算してplotする。
# Dmeanの計算は，別スクリプトであるcal_Dmean_by_sim.pyをimportして使っている。

import os
import argparse

import pandas as pd
import matplotlib #to set use('Agg') 
matplotlib.use('Agg') #place this before any other pylab/matplotlib/pyplot import
import matplotlib.pyplot as plt
import seaborn as sns

import cal_Dmean_by_sim

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--target_seqregion", type=str, help="target sequence region. sep=,")
parser.add_argument("-o", "--output_dir", type=str, default="../results/Dmean/pol-domains-HXB2-Pfam")
args = parser.parse_args()


def plot_barplot_Dmean(df, outfile='./out_plot_bar_Dmean.png'):
#    sns.set(style="whitegrid")
    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    fig, ax = plt.subplots(1, 1, figsize=(8,6)) #fig->figure obj. ax->graph obj. 2,1とかだとgは配列に.2,2だとarray.
    ax = sns.barplot(x='seq_region', y='Dmean', data=df)
    plt.xlabel('sequence region')
    plt.ylabel('Dmean')
    fig.savefig(outfile, dpi=180)
    sns.plt.close()

def plot_barplot_variance(df, outfile='./out_plot_bar_Variance.png'):
#    sns.set(style="whitegrid")
    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    fig, ax = plt.subplots(1, 1, figsize=(8,6)) #fig->figure obj. ax->graph obj. 2,1とかだとgは配列に.2,2だとarray.
    ax = sns.barplot(x='seq_region', y='Variance', data=df)
    plt.xlabel('sequence region')
    plt.ylabel('variance of distance')
    fig.savefig(outfile, dpi=180)
    sns.plt.close()


if __name__ == '__main__':
    if args.target_seqregion:
        seqregions = args.target_seqregion.split(',')
    else:
        seqregions = "RVP RVT_1 RVT_thumb RVT_connect RNase_H Integrase_Zn rve IN_DBD_C".split(' ')

    outdir = args.output_dir
    os.makedirs(outdir, exist_ok=True)
    out_file_Dmean = outdir + "/out_plot_bar_Dmean_pi90.png"
    out_file_Variance = outdir + "/out_plot_bar_Variance_pi90.png"


    in_dir = "../analysis/sim/pol-domains-HXB2-Pfam/"
    file_name_core = "_pi90_blastp_1e-5_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.txt"

    dmeans = []
    variances = []
    for seqre in seqregions:
        in_file = in_dir + seqre + file_name_core

        print(in_file)
        in_fh = open(in_file, 'r')
        dmean, variance = cal_Dmean_by_sim.cal_Dmean_var_file(in_fh)
        print("Dmean: ", dmean)
        print("Variance: ", variance)
        dmeans.append(dmean)
        variances.append(variance)

    df = pd.DataFrame({'seq_region':seqregions, 'Dmean':dmeans, 'Variance':variances})
    plot_barplot_Dmean(df, out_file_Dmean)
    plot_barplot_variance(df, out_file_Variance)
