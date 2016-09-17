#!/usr/bin/env python3

import argparse

import pandas as pd

# 2016.4.25 by Shohei Nagata
# データセットまとめ向け
# attrsファイルを読み込んで，fastaも読み込んで，使用した配列のみでデータセット表をdfで保持して，任意の情報を取り出せるようなscript


parser = argparse.ArgumentParser()
parser.add_argument("-a", "--input_attrs_file", type=str, help="Input Attributes file.", default="../analysis/attrs_HIV,SIV_regionyear.txt")
parser.add_argument("-f", "--input_fasta_file", type=str, help="Input FASTA file.", default="../data/A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta")
args = parser.parse_args()

def collect_fasta_ids(fasta_file):
    ids = []
    with open(fasta_file) as fasta_fh:
        for line in fasta_fh:
            if '>' in line:
                header = line.rstrip().split('>')[-1]
                ids.append(header)
    return ids

    



if __name__ == '__main__':
    df = pd.read_csv(args.input_attrs_file, delimiter='\t')
    ids = collect_fasta_ids(args.input_fasta_file)

    #FASTAで使われている配列のみ
    datasets_df = df[df['ID'].isin(ids)]
    #地道にやっていく.. 16.4.25
    for st in "A1 A2 B C".split():
        print("\n",st)
        print(datasets_df[df["Subtype"] == st].groupby("Region").size())
        print(datasets_df[df["Subtype"] == st].groupby("SamplingYear").size())
    quit()


    #FASTA関係なし 16.6.20
    datasets_df = df[df['ID'].isin(ids)]
    print("\n\nFASTA関係なしで")

    nums_HIV1_subtype = df[df["Organism"] == "HIV-1"].groupby("Subtype").size()
    nums_HIV1_subtype.to_csv("datasets_HIV-1_subtypes.csv")

    for st in "A1 A2 B C D F1 F2 G H J K U".split():
        print("\n",st)
        print(len(df[df["Subtype"] == st]))
        print(df[df["Subtype"] == st].groupby("Region").size())
        print(df[df["Subtype"] == st].groupby("SamplingYear").size())
    #

#ぼつ
#    print(datasets_df[["Subtype", "Region"]].groupby("Subtype").size())
#    print(datasets_df["Subtype"].value_counts())
#    print(datasets_df["Region"].value_counts())
#    print(datasets_df["SamplingYear"].value_counts())

