
#from math import log
import sys
from collections import Counter, defaultdict
import pprint

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import motifs
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.Alphabet import IUPAC, Gapped

import sn_utils

# 2016.8.17 memo
# 累積相対エントロピーのとこが未完，calc_plot_CRE_using_hmmbuild.py へ。

#2016.7.10 by Shohei Nagata.
#fastaファイルを読み込ませる。entoropyのplot機能と，KLダイバージェンスのplot機能，累積相対エントロピーのplot機能をのplot機能をつける
#ドメインなどの位置情報が書かれたtomlファイルも合わせて読み込ませると，その位置を表示してくれる

#サブタイプA,B,Cのみ対応なう。ハードコーディングだからだよ。

#sns.set('poster', 'whitegrid', 'dark', font_scale=1.2,
#        rc={"lines.linewidth": 2, 'grid.linestyle': '--'})
fig_width = 8
font_size = 20

print = pprint.pprint

#functionally same as Counter
def counter(seq):
    dic = {}
    for c in seq:
        try:
            dic[c] += 1
        except KeyError:
            dic[c] = 0
            dic[c] += 1
    return dic

def distribute(s):
#    print(list(counter(s).values()))
    return np.array(list(Counter(s).values()))

def distribute_fitshape(s1,s2):
    s1c_d = Counter(s1)
    s2c_d = Counter(s2)

    return np.array(list(Counter(s).values()))

#未修正，
#def ori_entropy(s): #scipy使うより自前で計算した方が計算結果桁数多い
#    p = counter(s)
#    lns = float(len(s))
#    return -sum( count/lns * log(count/lns,2) for count in p.values())
#

def calc_entropy(p):
    return stats.entropy(p, base=2)  

def calc_kl_divergence(p, q):
    return stats.entropy(p, q, base=2)  

def plot_nparray(array, outfile="./out_entropy.png", ylabel=""):
    fig, ax = plt.subplots(1,1, figsize=(fig_width, fig_width*0.6))
    ax.plot(array, c='k')
    ax.set_ylabel('Entropy', fontsize=font_size)
    fig.savefig(outfile, dpi=180)
    sns.plt.close()

def plot_nparray_with_annopos(array, annopos_d, outfile="./out_entropy_annopos.png", ylabel=""):
    sns.set('poster', 'whitegrid', 'dark', font_scale=1.2,
            rc={"lines.linewidth": 2, 'grid.linestyle': '--'})

    fig, axs = plt.subplots(2,1, 
                            sharex=True, 
                            figsize=(fig_width, fig_width*0.6),
                            gridspec_kw={'height_ratios':[6, 1]})

    #plot entropy
    ax = axs[0]
    ax.plot(array, c='k')
    ax.set_ylabel('Entropy', fontsize=font_size)

    ax = axs[1]
    sns.set_style('white')
    sn_utils.draw_annopos(ax, annopos_d)
    ax.set_yticks([])
    ax.set_xlabel('Position', fontsize=font_size)

    # Final touches
    plt.tight_layout(rect=(0.0, 0.02, 0.98, 0.98), pad=0.1, h_pad=0.5, w_pad=0.4)
    fig.savefig(outfile, dpi=180)
    sns.plt.close()

def smooth_array(array, window_size, min_valid_fraction=0.30):
    smoothed_divergence = sn_utils.running_average_masked(array, window_size, min_valid_fraction=min_valid_fraction)
    return smoothed_divergence
    
def fasta_entropy_a(fasta_text):
    fasta_d = sn_utils.read_fasta(fasta_text)
    seqs = []
    for v in fasta_d.values():
        seqs.append(v.seq)

    seqlen = len(seqs[0])
    entropy_l = []
    for i in range(seqlen):
        col = ""
        for seq in seqs:
            col += seq[i]
        entropy_num = calc_entropy(distribute(col))
        entropy_l.append(entropy_num)
#        print(i,entropy_num)
    return np.array(entropy_l)

#def fasta_cre_a_bio(fasta_text):  #未完成
#    """
#    Cumulative Relative Entropy
#    # Biopythonの利用を試みたが失敗。ただの残骸コードなう
#    # 一旦fastaをdictへしてから配列取り出してまたalignment objectにしているのが無駄無駄
#    """
#    alphabet_protein = Gapped(IUPAC.protein)
#
#    fasta_d = sn_utils.read_fasta(fasta_text)
#    grouped_fasta_d = group_fasta_seqs(fasta_d)
##    pprint.pprint(grouped_fasta_d);quit()
#    seqs_d = defaultdict(list)
#    cre_a = np.array([])
#
#    #以下冗長コード
#    for t, fa_d in grouped_fasta_d.items():
#        for v in fa_d.values():
#            seqs_d[t].append(v.seq) #SeqRecord としているので注意。以前は.seqで配列のみ取り出していた。
#
#    for t_main, seqs_main in seqs_d.items():
#        seqlen = len(seqs_main[0])
#        entropy_l = []
#        seqs_others = []
#        for t_others, seqs_tmp in seqs_d.items():
#            if t_others != t_main:
#                seqs_others.append(seqs_tmp)
#
#        print(alphabet_protein)
#        m_main = motifs.create(seqs_main, alphabet=alphabet_protein)
#        print(m_main);quit()
#        
#        summary_align_main = AlignInfo.SummaryInfo( MultipleSeqAlignment(seqs_main) )
#        print(summary_align_main);quit()
#        summary_align_others = AlignInfo.SummaryInfo( MultipleSeqAlignment(seqs_others) )
#
#        for i in range(seqlen):
#            col_main = ""
#            col_others = ""
#            for seq in seqs_main:
#                col_main += seq[i]
#            for seq in seqs_others:
#                col_others += seq[i]
#
#            print(distribute(col_main))
#            print(distribute(col_others))
#            entropy_num = calc_kl_divergence(distribute(col_main), distribute(col_others))
#            entropy_l.append(entropy_num)
#    #        print(i,entropy_num)
#        cre_a = cre_a + np.array(entropy_l)
#        pprint.pprint(cre_a); quit()
#
#    return cre_a
#    
#def fasta_cre_a(fasta_text):   #動きまへん。
#    """
#    Cumulative Relative Entropy
#    """
#    fasta_d = sn_utils.read_fasta(fasta_text)
#    grouped_fasta_d = group_fasta_seqs(fasta_d)
##    pprint.pprint(grouped_fasta_d);quit()
#    seqs_d = defaultdict(list)
#    cre_a = np.array([])
#
#    #以下冗長コード
#    for t, fa_d in grouped_fasta_d.items():
#        for v in fa_d.values():
#            seqs_d[t].append(v.seq)
#
#    for t_main, seqs_main in seqs_d.items():
#        seqlen = len(seqs_main[0])
#        entropy_l = []
#        seqs_others = []
#        for t_others, seqs_tmp in seqs_d.items():
#            if t_others != t_main:
#                seqs_others.append(seqs_tmp)
#
#        for i in range(seqlen):
#            col_main = ""
#            col_others = ""
#            for seq in seqs_main:
#                col_main += seq[i]
#            for seq in seqs_others:
#                col_others += seq[i]
#
#            print(distribute(col_main))
#            print(distribute(col_others))
#            entropy_num = calc_kl_divergence(distribute(col_main), distribute(col_others))
#            entropy_l.append(entropy_num)
#    #        print(i,entropy_num)
#        cre_a = cre_a + np.array(entropy_l)
#        pprint.pprint(cre_a); quit()
#
#    return cre_a
#
#def group_fasta_seqs(fasta_d):
#    """
#    サブタイプごとにわける
#    ハードコーディングなう.もっと良い方法求む
#    """
#    grouped_d = defaultdict(dict)
#    target_list = [":A",":B",":C"]
#    for k, v in fasta_d.items():
#        for t in target_list:
#            if t in k:
#                grouped_d[t].update({k:v})
#    return grouped_d


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--input_toml_file", type=str, help="Input TOML file.", default="/Users/NagataShohei/Documents/bio-study/16spp/scripts/toml/aligned_HIV-1_gag-pol_Pfam_HXB2.toml")
    parser.add_argument("-f", "--input_fasta_file", type=str, help="Input FASTA file.")
    parser.add_argument("-s", "--seq_type", type=str, default="aa")
    args = parser.parse_args()

    seq_type = args.seq_type
    if seq_type == 'nuc':
        min_valid_fraction = 0.95
    elif seq_type == 'aa':
        # Two out of three are masked by design
        min_valid_fraction = 0.30

    
    annopos_d = sn_utils.read_tomlfile(args.input_toml_file)

#    #cre
#    cre_ma = np.ma.array(fasta_cre_a(sys.stdin))
#    quit()



    #entropy
    entropy_ma = np.ma.array(fasta_entropy_a(sys.stdin))
    s_entropy_ma = smooth_array(entropy_ma, 50)

    plot_nparray_with_annopos(s_entropy_ma, annopos_d, outfile="./out_entropy_annopos_mv50.png")
    plot_nparray_with_annopos(entropy_ma, annopos_d, outfile="./out_entropy_annopos.png")

#    plot_nparray(s_entropy_ma)


        

        

            

