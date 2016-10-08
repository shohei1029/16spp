#!/usr/bin/env python3

import time
import argparse

from progressbar import ProgressBar

#memo
#for all-against-all BLAST
#->全てのqueryがどっかしらにhitする前提
#all-against-all以外でも使えるようにするには，queryとして使われた配列がinput fastaの何番目の配列かってのをみないといけない
#  全部のqueryがどっかにhitするわけじゃないから。

#毎回出力ファイルを最初から全部読むんじゃなくて，tail -f とか使って効率化していきたい

parser = argparse.ArgumentParser()
parser.add_argument("-i","--blast_inputfile",type=str,help="Blast input file.")
parser.add_argument("-o","--blast_outputfile",type=str,help="Blast output file.")
args = parser.parse_args()


with open(args.blast_inputfile, 'r') as bl_in_fh:
    lines = bl_in_fh.read()
    numentry = lines.count('>')

pbar = ProgressBar(max_value=numentry,redirect_stderr=True)
pbar.start()

while True:
    with open(args.blast_outputfile, 'r') as bl_out_fh: 
        lines = bl_out_fh.read()
        prog_numquery = lines.count('# Query: ')

        pbar.update(prog_numquery)
        time.sleep(3)

        if prog_numquery >= numentry:
            break

pbar.finish()
        

