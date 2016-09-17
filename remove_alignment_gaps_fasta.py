#!/usr/bin/env python3

import sys

from Bio import SeqIO

import sn_utils

# 2016.9.1 by Shohei Nagata
# Remove from an alignment all columns that contain a gap in the specified species.
# Usage: remove_alignment_gaps_fasta.py <alignment> <species; fasta header>

# headerのスペース以前のみでリファレンスを指定しているので，スペース以前の文字列が同じというheaderが複数含まれていたらやばい。

def remove_char_at_index(string, index):
    new_string =  string[:index] + string[index+1:]
    return new_string

if __name__ == '__main__':
    fasta_d = sn_utils.read_fastafile(sys.argv[1]) #fasta headerのスペース以降は消える。
    ref_id  = sys.argv[2].split(' ')[0] #スペース前のみに合わせる
    ref_seq = fasta_d[ref_id].seq
#    print(ref_seq)

    gap_pos = []
    for i,c in enumerate(ref_seq):
        if c == '-':
            gap_pos.append(i)

    #n^2 ..orz
    for h, seqrcd in fasta_d.items():
        for j,p in enumerate(gap_pos):
            p = p - j 
            seqrcd.seq = remove_char_at_index(seqrcd.seq, p)
        SeqIO.write(seqrcd, sys.stdout, 'fasta')
