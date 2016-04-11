#!/usr/bin/env python3

import sys
import time
import os
import argparse

import matplotlib #to set use('Agg') 
matplotlib.use('Agg') #place this before any other pylab/matplotlib/pyplot import
import matplotlib.pyplot as plt
import seaborn as sns


# 2016.3.7
# Shohei N.
#
# sys.stdinになんかのリスト流し込むとヒストグラムにして出力してくれる。おわり。

parser = argparse.ArgumentParser()
parser.add_argument("-b","--bin", type=int, help="bin of histgram (int)")
parser.add_argument("-d","--kde", type=bool, default=False, help="Whether to plot a gaussian kernel density estimate.")
args = parser.parse_args()

def plot_distplot(x, outfile="./out_plot_hist.png", b=args.bin, kde=args.kde):
#    sns.set(style="whitegrid")
    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    fig, ax = plt.subplots(1, 1, figsize=(8,6)) #fig->figure obj. ax->graph obj. 2,1とかだとgは配列に.2,2だとarray.
    ax = sns.distplot(x, bins=b, kde=kde)
    fig.savefig(outfile, dpi=180)
    sns.plt.close()
    

if __name__ == '__main__':
    l = []
    for line in sys.stdin:
        line = line.rstrip()
        l.append(float(line))
    plot_distplot(l)
