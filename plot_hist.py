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
# sys.stdinになんかのリスト流し込むとヒストグラムにして出力してくれる。おわり。

# 2016.5.26
# x軸の範囲とか，タイトルとか，引数で全部与えられるように。

parser = argparse.ArgumentParser()
parser.add_argument("-b","--bin", type=int, help="bin of histgram (int)")
parser.add_argument("-d","--kde", type=bool, default=False, help="Whether to plot a gaussian kernel density estimate.")
parser.add_argument("-o","--output_file", type=str)
parser.add_argument("--xlim", type=str, help="sep=,   ex. 0,1")
parser.add_argument("--ylim", type=str, help="sep=,   ex. 0,10000")
parser.add_argument("--xlabel", type=str)
parser.add_argument("--ylabel", type=str)
parser.add_argument("--title", type=str)
args = parser.parse_args()

def plot_distplot(x, outfile=args.output_file, b=args.bin, kde=args.kde, xlabel=args.xlabel, ylabel=args.ylabel, xlim_str=args.xlim, ylim_str=args.ylim, title=args.title):
#    sns.set(style="whitegrid"k
    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    fig, ax = plt.subplots(1, 1, figsize=(8,6)) #fig->figure obj. ax->graph obj. 2,1とかだとgは配列に.2,2だとarray.
    ax = sns.distplot(x, bins=b, kde=kde)

    if xlim_str:
        xlim_l = xlim_str.split(",")
        ax.set_xlim(int(xlim_l[0]), int(xlim_l[1]))
    if ylim_str:
        ylim_l = ylim_str.split(",")
        ax.set_ylim(int(ylim_l[0]), int(ylim_l[1]))

    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if title:
        ax.set(title=title)

    if not outfile:
        outfile = "./out_plot_hist.png"

    fig.savefig(outfile, dpi=300)
    sns.plt.close()
    print("saved as", outfile)
    

if __name__ == '__main__':
    l = []
    for line in sys.stdin:
        line = line.rstrip()
        l.append(float(line))

    plot_distplot(l)
