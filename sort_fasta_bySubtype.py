
import sys
from collections import defaultdict

from Bio import SeqIO

import sn_utils

# 2016.8.4 by Shohei Nagata
# sort fasta sequences by subtype (in fasta header) to use "multi-Harmony: multi-group Sequence Harmony & multi-Relief"
# fasta -> stdin, sorted fasta -> stdout, a list of group sizes -> stderr


#計算量無駄ェ..
if __name__ == "__main__":
    fasta_d = sn_utils.read_fasta(sys.stdin)
    subtypes_d = defaultdict(dict)
    subtypes_s = set()

    for k, v in fasta_d.items():
        subtype = k.split(":")[-1]
        subtypes_d[subtype].update({k:v})

    gs = ""
    for s, seqd in sorted(subtypes_d.items()):
        for h, seqobj in seqd.items():
            SeqIO.write(seqobj, sys.stdout, "fasta")
        gs += str(len(seqd)) + " "

    print(gs, file=sys.stderr)
