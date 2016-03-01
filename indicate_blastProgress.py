#!/usr/bin/env python3

import time
import argparse

from progressbar import ProgressBar

parser = argparse.ArgumentParser()
parser.add_argument("-i","--blast_inputfile",type=str,help="Blast input file.")
parser.add_argument("-o","--blast_outputfile",type=str,help="Blast output file.")
args = parser.parse_args()


with open(args.blast_inputfile, 'r') as bl_in_fh:
    lines = bl_in_fh.read()
    numentry = lines.count('>')

pbar = ProgressBar(max_value=numentry,redirect_stderr=True)
pbar.start()

#先頭に#かあるがみて，そのなかでカウントする
while True:
    with open(args.blast_outputfile, 'r') as bl_out_fh:
        print("with")
        lines = bl_in_fh.read()
        prog_numquery = lines.count('# Query: ')
        print(prog_numquery)

        pbar.update(prog_numquery)

    if prog_numquery >= numentry:
        print("break")
        break


pbar.finish()
        

