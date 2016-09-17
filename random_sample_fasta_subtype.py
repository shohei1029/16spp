
import sys
from collections import defaultdict

from Bio import SeqIO

import sn_utils

# 2016.9.2 by Shohei Nagata
# sys.stdinでfastaを読み込ませ，各サブタイプごとにランダムに指定数配列とってきて出力する
# サブタイプ情報はheaderをもとに。ドメイン領域(subtype:A-hoge)みたいに書いてあっても対応している

seqnum = 10

if __name__ == '__main__':
    subtype_d = sn_utils.read_fasta_subtype(sys.stdin)

    # merge subtype A1 and A2, F1 and F2 (hard coding..)
    if 'A2' in subtype_d.keys():
        subtype_d['A'] = {}
        subtype_d['A'].update(subtype_d['A1'])
        subtype_d['A'].update(subtype_d['A2'])
        del subtype_d['A1']
        del subtype_d['A2']
    if 'F2' in subtype_d.keys():
        subtype_d['F'] = {}
        subtype_d['F'].update(subtype_d['F1'])
        subtype_d['F'].update(subtype_d['F2'])
        del subtype_d['F1']
        del subtype_d['F2']

    sampled_fasta_d = defaultdict(dict)
    for k,v in sorted(subtype_d.items()):
        sampled_fasta_d[k] = sn_utils.random_sample_fasta(v, seqnum)

        for seqobj in sampled_fasta_d[k].values():
            sys.stdout.write(">{}\n{}\n".format(seqobj.id, seqobj.seq))
    

    
