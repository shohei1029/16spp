#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict

from Bio import SeqIO
import toml

# 2016.5.12, Shohei Nagata
# swiss-protを基にしたtomlに書かれた位置を，アライメント後fastaの1番目の配列上の位置へ変換して(ギャップ分変換)，新しくtomlを吐き出すスクリプト.
# 出力は標準出力

#memo

parser = argparse.ArgumentParser()
parser.add_argument("-t","--input_toml_file", type=str, help="Input TOML file.")
parser.add_argument("-f","--input_fasta_file", type=str, help="Input FASTA file.")
args = parser.parse_args()


def seqpos_to_gapped_seqpos(pos,seq):
    c_count = 0
    gapped_pos = 0
    for c in seq:
        if c != '-':
            c_count += 1
        gapped_pos += 1
        if c_count == pos:
            return gapped_pos

def read_1st_seq_fasta(fasta_file):
    with open(fasta_file,'r') as fa_fh:
        for seqrcd in SeqIO.parse(fa_fh,"fasta"):
            return seqrcd.seq

def read_toml(toml_file):
    with open(toml_file) as tomlfh:
        toml_d = toml.loads(tomlfh.read())
    return toml_d

def convert_toml(toml_d,seq):
    new_seqregion_d = dict()
    for k,l in toml_d.items():
        new_start_pos = seqpos_to_gapped_seqpos(l[0],seq)
        new_end_pos   = seqpos_to_gapped_seqpos(l[1],seq)
        new_seqregion_d.update({k: [new_start_pos,new_end_pos]})
    return new_seqregion_d

def write_toml(toml_d):
    sys.stdout.write(toml.dumps(toml_d))



if __name__ == '__main__':
    toml_file = args.input_toml_file
    in_fasta_file = args.input_fasta_file

    in_toml_d = read_toml(toml_file)
    stdseq = read_1st_seq_fasta(in_fasta_file)

    seqregion_d = convert_toml(in_toml_d,stdseq)
    write_toml(seqregion_d)
