#!/usr/bin/env python

import sys
import random

from Bio import SeqIO

#2016.4.1 by S.N.
#fastaファイルを読み込ませ，指定した数のみランダムサンプリングして書き出す
#ちゃんとrandom関数を使用している(辞書では順序が保存されないってのがランダムなのかわからないため)

n = 1000

def read_fasta_todict(fasta_txt):
    fa_dict = SeqIO.to_dict(SeqIO.parse(fasta_txt, "fasta"))
    return fa_dict

def random_sample_fasta(fa_d, n):
    fa_sampled_l = random.sample(list(fa_d.items()), n)
    fa_sampled_d = dict()
    for seq in fa_sampled_l: #データ形式をSeqIOで読み込んだものと同じに戻す
        fa_sampled_d.update({seq[0]: seq[1]})
    return fa_sampled_d

def write_fasta_seqio_list(fa_d, out_fh=sys.stdout): 
    for seqobj in fa_d.values():
        out_fh.write(">{}\n{}\n".format(seqobj.id, seqobj.seq))

if __name__ == '__main__':
    fa_d = read_fasta_todict(sys.stdin)
    fa2_d = random_sample_fasta(fa_d, n)
    write_fasta_seqio_list(fa2_d)
