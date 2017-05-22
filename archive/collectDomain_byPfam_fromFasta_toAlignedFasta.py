

import os
import sys
import time
import re
from collections import defaultdict
from Bio import SeqIO
import multiprocessing as mp
#import subprocess

stime = time.time()
argvs = sys.argv

#process number
#proc = 8
proc = int(argvs[2]) if argvs[2] else 4

#target_domains_name :Pfam

#############
# I/O files #
#############
hmmscan_file = "./Pfam_HMMER/Pfam-hmmscan_HIV-1-gM-noRs_pol-aa_v3.txt" 

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


def find_all_files(directory):
    for root, dirs, files in os.walk(directory):
        #yield root
        for file in files:
            yield os.path.join(root, file)



#####################
# Read hmmscan file #
#####################
domain_dict = defaultdict(lambda:defaultdict(lambda:defaultdict()))
#data_dict = defaultdict(lambda:defaultdict(lambda:defaultdict()))

domain_name_dict = defaultdict(str) #set的な用途で。
with open(hmmscan_file,'r') as hmmscan_fh:
    for line in hmmscan_fh:
        if line[0] == '#':
            continue
        line = line.rstrip()

        p = re.compile('\s+')
        hmmscandata = p.split(line)

        domain_name  = hmmscandata[0]
        #query_acc    = hmmscandata[3].split('|')[0]
        query_acc    = hmmscandata[3]
        domain_start = int(hmmscandata[17])
        domain_end   = int(hmmscandata[18])
        domain_des   = (' ').join(hmmscandata[22:])

        domain_name_dict[domain_name] = domain_des
        domain_dict[query_acc][domain_name]["start"] = domain_start
        domain_dict[query_acc][domain_name]["end"]   = domain_end

#print(domain_dict)

#print Domains
#for dk,dv in domain_name_dict.items():
#    print(dk,": ",dv)


##################################################################
# Search and read MAFFT aligned FASTA files (gap included fasta) #
##################################################################
# search
target_files = []
targetdir = argvs[1]
for file in find_all_files(targetdir):
    if "mafft" in file.split('/')[-1] : #read MAFFT aligned files only
        ext = file.split('/')[-1].split('.')[-1] 
        if ext=="fasta" or ext=="faa" or ext=="fna" or ext=="fa":
            #target_files.add(file)
            target_files.append(file)


# read and process
#並列化できるぜ！
#for fasta_file in target_files:
def make_fasta_by_domain(fasta_file):
    print("processing..",fasta_file)
    acc_set = set()
    #emblacc_set = set()
    with open(fasta_file,'r') as fa_fh:
        fa_dict = SeqIO.to_dict(SeqIO.parse(fa_fh, "fasta"))

#    fa2_dict = defaultdict(None) #キーがaccのみ。(もとはacc|HIV-1|subype:hoge みたいになってた)
#    for k,v in fa_dict.items():
#        tmp_acc = k.split('|')[0]
#        acc_set.add(tmp_acc)
#        fa2_dict[tmp_acc] = v #accのみをキーとしたdictへ。
    fa2_dict = fa_dict #キーをもとのやつに戻すｗ

    fastafile_name = fasta_file.split('/')[-1].split('.')[0]
    domain_outdir = "/".join(fasta_file.split('/')[0:-1]) + "/domain/" + fastafile_name
#    print(domain_outdir)
    os.makedirs(domain_outdir,exist_ok=True) 



#strするときにpos-1で指定することを忘れずに！
#まずドメインごとのsequence logo 出力を実装する


    for domn in domain_name_dict.keys():
        domain_fasta = "{}/{}.fasta".format(domain_outdir,domn)
        dom_start,dom_end = 100000000000,0
        #print(domain_fasta)

#        with open(domain_fasta,'w') as dofa_fh:
#            #先に最長となるようなDomainの開始・終了点を決めて，次のfor文で配列とってきて書き込んでいく
#            for k_acc in fa2_dict.keys():
#                try:
#                    dom_start_tmp = seqpos_to_gapped_seqpos(domain_dict[k_acc][domn]["start"],fa2_dict[k_acc].seq)
#                    dom_end_tmp   = seqpos_to_gapped_seqpos(domain_dict[k_acc][domn]["end"],fa2_dict[k_acc].seq)
#
#                    dom_start = dom_start_tmp if dom_start_tmp < dom_start else dom_start
#                    dom_end   = dom_end_tmp   if dom_end_tmp   > dom_end   else dom_end
#                except KeyError as err: #主にそのドメインがヒットしなかった配列
#                    pass
#            for k_acc in fa2_dict.keys():
#                fasta_header = "{}-{}".format(k_acc,domn)
#                try:
#                    dom_seq   = fa2_dict[k_acc].seq[dom_start-1:dom_end-1]
#                except KeyError as err: #主にそのドメインがヒットしなかった配列
#                    pass
#                dofa_fh.write(">{}\n{}\n".format(fasta_header,dom_seq))

        #先に最長となるようなDomainの開始・終了点を決めて，次のfor文で配列とってきて書き込んでいく
        for k_acc in fa2_dict.keys():
            try:
                dom_start_tmp = seqpos_to_gapped_seqpos(domain_dict[k_acc][domn]["start"],fa2_dict[k_acc].seq)
                dom_end_tmp   = seqpos_to_gapped_seqpos(domain_dict[k_acc][domn]["end"],fa2_dict[k_acc].seq)

                dom_start = dom_start_tmp if dom_start_tmp < dom_start else dom_start
                dom_end   = dom_end_tmp   if dom_end_tmp   > dom_end   else dom_end
            except KeyError as err: #主にそのドメインがヒットしなかった配列
                pass
        print("{}\t{}\t{}".format(domn,dom_start,dom_end))



##########
# main() #
##########
#Multiprocessing
pool = mp.Pool(proc)
cb = pool.map(make_fasta_by_domain,target_files)

print(time.time() - stime,"[s]")
