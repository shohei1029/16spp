#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO


# v2 (16spp)
#change log
# 2016.2.29
#  copied from 15fall v1
#  stdin -> stdout (default)

#memo
#nt 未実装

class END(Exception):
    pass

# Argument
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",type=str,help="Input fasta file.")
parser.add_argument("-o","--output_file",type=str,help="Output fasta file.")
parser.add_argument("-t","--seq_type",type=str,choices=["aa","nt"],default="aa",help="Fasta sequenece type.") #nt 未実装
args = parser.parse_args()

# prep
aa_x_set = {"B","Z","X","J"}

if __name__ == '__main__':
    if args.input_file:
        in_fh = open(args.input_file,"r")
    else:
        in_fh = sys.stdin
    
    if args.output_file:
        out_fh = open(args.output_file,"w")
    else:
        out_fh = sys.stdout
    

    fastaobj = SeqIO.to_dict(SeqIO.parse(in_fh,"fasta"))
    
    count = 0
    for k,v in fastaobj.items():
        try:
            for x in aa_x_set:
                if x in set(v.seq):
                    count += 1
                    raise END
            out_fh.write(">{h}\n{s}\n".format(h=v.id,s=v.seq))
        except END:
            pass
    
    sys.stderr.write(" ".join(["Delete",str(count),"seqs (",len(fastaobj),")"]))
    sys.stderr.write(" ".join(["Delete",str(count/len(fastaobj)*100),"%"]))
