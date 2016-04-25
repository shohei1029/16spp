#!/usr/bin/env python3

import sys
import argparse
import statistics

from Bio import SeqIO

# v1 (16spp)
# change log
# 2016.4.9
#  created file based on filter_fasta_seq.py
#  stdin -> stdout (default)
# 2016.4.13
#  fastaをdictとして読み込む時の，ValueError: Duplicate key 'DQ164129|HIV-1|subtype:C|gag'への対処(例外処理)
#  ->同じIDがあるときは2個め無視！その分配列数減るので注意！
#(まだしてない↑)

#memo
#基準とする配列長の指定%未満の配列を取り除いて出力する。
#-n optionでfasta headerの最後にある遺伝子名/タンパク質名のとこも除去する。(./getseq_ingbk_withbginfo_tofasta.py で出力すると遺伝子名/タンパク質名が最後につくため)

#class END(Exception):
#    pass

# Argument
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file", type=str, help="Input fasta file.")
parser.add_argument("-o","--output_file", type=str, help="Output fasta file.")
parser.add_argument("-r","--ref_seq_length", type=int, default=1000, help="Reference sequence length.") #nt 未実装
parser.add_argument("-c","--criteria", type=float, default=0.9, help="remove sequence length under (critetia*100)% of reference sequence length.") 
parser.add_argument("-n","--remove_gene_name", type=bool, default=False, help="remove gene name in fasta headers.") 
args = parser.parse_args()


if __name__ == '__main__':
    if args.input_file:
        in_fh = open(args.input_file,"r")
    else:
        in_fh = sys.stdin
    
    if args.output_file:
        out_fh = open(args.output_file,"w")
    else:
        out_fh = sys.stdout
    

    try:
        fastaobj = SeqIO.to_dict(SeqIO.parse(in_fh,"fasta"))
    except ValueError:
        pass #これでDQ164129にgagが2つあることで発生するValueError: Duplicate key 'DQ164129|HIV-1|subtype:C|gag'を防げるのかはわからない。2こめのgagは無視してよい。
        
    
    seqlens_all = set()
    seqlens_filterd = set()
    cri_len = args.ref_seq_length * args.criteria
    for k,v in fastaobj.items():
        slen = len(v.seq)
        #print(slen)
        seqlens_all.add(slen)
        if slen < cri_len:
            continue
        seqlens_filterd.add(slen)

        if args.remove_gene_name:
            fah = "|".join(v.id.split("|")[:-1])
        else:
            fah = v.id

        out_fh.write(">{h}\n{s}\n".format(h=fah,s=v.seq))

    
    #puts infomation
    sys.stderr.write("Filtering fasta by sequence length..\n")
    sys.stderr.write("Before: \n")
    sys.stderr.write("mean: " + str(statistics.mean(seqlens_all)) + "\n")
    sys.stderr.write("median: " + str(statistics.median(seqlens_all)) + "\n")
    sys.stderr.write("After: \n")
    sys.stderr.write("mean: " + str(statistics.mean(seqlens_filterd)) + "\n")
    sys.stderr.write("median: " + str(statistics.median(seqlens_filterd)) + "\n")
    
