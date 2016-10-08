#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict

from Bio import SeqIO
import toml

# 2016.4.24 by Shohei Nagata
# そこのドメインの位置情報をもとに，配列を切り出し，ドメインごとのfastaへしてstdoutへするスクリプト。
# ドメインの位置情報をPfamからとってきてtomlに書いたものを読みこませる。
# アライメント済のfastaファイルを読みこませる。ちゃんとgapを抜いて配列を切り出してくれる.
# 引数にtomlファイル

# 2016.5.12 メモ
# ..これって，アライメントしたときHXB2配列にギャップが入ってない前提になってないか？
# ☆このスクリプトはそのままで，toml側をHXB2のギャップあり位置へ変換するスクリプトを作成する.

#memo

parser = argparse.ArgumentParser()
parser.add_argument("-t","--input_toml_file", type=str, help="Input TOML file.")
parser.add_argument("-f","--input_fasta_file", type=str, help="Input FASTA file.")
args = parser.parse_args()


#############
# I/O files #
#############

#############
# Functions #
#############
def seqpos_to_gapped_seqpos(pos,seq):
#    '''
#    ギャップなしFASTAファイルのpositionをギャップ有りFASTAファイルのpositionへ変換
#    convert char position in ungapped FASTA to gapped FASTA. (I'm sorry I can't use English well.)
#    '''
    c_count = 0
    gapped_pos = 0
    for c in seq:
        if c != '-':
            c_count += 1
        gapped_pos += 1
        if c_count == pos:
            return gapped_pos


def read_fasta_todict(fasta_file):
    with open(fasta_file,'r') as fa_fh:
        fa_dict = SeqIO.to_dict(SeqIO.parse(fa_fh, "fasta"))
    return fa_dict

        
def make_seqsec_fasta(seq_secs_ld, alnfasta_dict):
    fa2_dict = alnfasta_dict

    for k,v in fa2_dict.items(): 
        sys.stderr.write("sequence length with gap: {}\n".format(len(v.seq)))
        break

#strするときにpos-1で指定することを忘れずに！
    for seqsec_name, spos_l in seq_secs_ld.items():
        for k_acc in fa2_dict.keys():
            fasta_header = "{}-{}".format(k_acc, seqsec_name)
            sec_seq   = str(fa2_dict[k_acc].seq[spos_l[0]-1:spos_l[1]-1])
            sec_seq_nogap = sec_seq.replace('-','')
            outfa_fh.write(">{}\n{}\n".format(fasta_header,sec_seq_nogap))


def read_seq_sec_fromToml(toml_file):
    with open(toml_file) as tomlfh:
        toml_d = toml.loads(tomlfh.read())
    return toml_d


##########
# main() #
##########
if __name__ == '__main__':
    toml_file = args.input_toml_file
    in_fasta_file = args.input_fasta_file
    outfa_fh = sys.stdout

    seq_sec = read_seq_sec_fromToml(toml_file)
    fa_d = read_fasta_todict(in_fasta_file)
    make_seqsec_fasta(seq_sec, fa_d)
