#!/usr/bin/env python3

import sys

import pandas as pd
import matplotlib #to set use('Agg') 
matplotlib.use('Agg') #place this before any other pylab/matplotlib/pyplot import
import matplotlib.pyplot as plt
import seaborn as sns

# Created on 2016.12.14 by Shohei Nagata.
# クラスターとサブタイプの一致率csv用に調整.
#  描写パラメータの調整。figsizeも，変えてしまうと文字が勝手に回転してから表示されるようになったりするよ。
# concordance_cluster_subtype_hoge.csv

def plot_heatmap(df, outfile=None, figsize=(8,16), annot=False, fmt='.2g', linewidths=0): #fmt:有効数字
#    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    fig, ax = plt.subplots(1, 1, figsize=figsize) #fig->figure obj. ax->graph obj. 2,1とかだとgは配列に.2,2だとarray.
    ax = sns.heatmap(df, cmap="YlGnBu", annot=annot, fmt=fmt, linewidths=linewidths)

    if not outfile:
        outfile = "./out_plot_heatmap.png"

    fig.savefig(outfile, dpi=180)
    sns.plt.close()
    print("saved as", outfile)


if __name__ == "__main__":
    df = pd.read_csv(sys.stdin)
    df = df.drop("seq identity threshold", axis=1)
#    plot_heatmap(df, figsize=(8,16), outfile="../analysis/heatmap_concordance_cluster_subtype.png")
    plot_heatmap(df, figsize=(8,19), annot=True, linewidths=0.2, outfile="../analysis/heatmap_concordance_cluster_subtype_annot.png")

#    df_t = df.T
#    plot_heatmap(df_t, figsize=(20,8), outfile="../analysis/heatmap_concordance_cluster_subtype_t.png")
