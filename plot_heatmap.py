#!/usr/bin/env python3

import sys

import pandas as pd
import matplotlib #to set use('Agg') 
matplotlib.use('Agg') #place this before any other pylab/matplotlib/pyplot import
import matplotlib.pyplot as plt
import seaborn as sns

# Created on 2016.12.14 by Shohei Nagata.
# 標準入力でcsvファイルを入力しheatmapを出力する.

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
    plot_heatmap(df)

