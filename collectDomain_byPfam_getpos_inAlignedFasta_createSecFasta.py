#!/usr/bin/env python3

import os
import sys
import time
import re
from collections import defaultdict
from Bio import SeqIO
#import multiprocessing as mp
#import subprocess

#memo
#Hmmerの結果ファイルからドメインの場所取得する
#アライメント後fataファイルを読み込ませ，ドメインの位置を数字で表示（最長に合わせる）
#手動で,ドメインもとにpolを４分割したときの境目数字を計算し，ハードコーディングする
#分割したときの各セクションについて，始め：終わりのリストの辞書を作成しておき，それを読み込ませて特定領域だけ切り出す。
# →ドメインやる場合も同じ方式でできるように！！すばらしい.←ドメインの配列長をドメインごとに揃えるとき限定な

# 2016.3.7 ドメインの配列長をドメインごとに揃えないようにする。そういう関数も作った。

stime = time.time()
#argvs = sys.argv

#process number
#proc = 8
#proc = int(argvs[2]) if argvs[2] else 4

#target_domains_name :Pfam

#############
# I/O files #
#############
hmmscan_file = "../data/Pfam-hmmscan_HIV-1-gM-noRs_pol-aa_v3.txt" 
in_alnfasta = "../data/mafft-linsi_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta"

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


domain_dict = defaultdict(lambda:defaultdict(lambda:defaultdict()))
domain_name_dict = defaultdict(str) #set的な用途で。
def read_hmmscan_domtblout_todict(hmmscan_file): #void func
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


def read_fasta_todict(fasta_file):
    with open(fasta_file,'r') as fa_fh:
        fa_dict = SeqIO.to_dict(SeqIO.parse(fa_fh, "fasta"))
    return fa_dict

#ドメイン領域の場所が入った辞書を返す。ドメイン領域は最長となるようにし，長さ揃える。
def get_domain_gpdpos_widest(alnfasta_dict):
    domain_gpdpos_ld = {}
    for domn in domain_name_dict.keys():
        dom_start,dom_end = 100000000000,0

        #先に最長となるようなDomainの開始・終了点を決めて，次のfor文で配列とってきて書き込んでいく
        for k_acc in alnfasta_dict.keys():
            try:
                dom_start_tmp = seqpos_to_gapped_seqpos(domain_dict[k_acc][domn]["start"],alnfasta_dict[k_acc].seq)
                dom_end_tmp   = seqpos_to_gapped_seqpos(domain_dict[k_acc][domn]["end"],alnfasta_dict[k_acc].seq)
                print(dom_end_tmp) #test

                dom_start = dom_start_tmp if dom_start_tmp < dom_start else dom_start
                dom_end   = dom_end_tmp   if dom_end_tmp   > dom_end   else dom_end
            except KeyError as err: #主にそのドメインがヒットしなかった配列
                pass
        domain_gpdpos_ld[domn] = [dom_start, dom_end]
    return domain_gpdpos_ld
    
#最長となるようなドメイン領域をとってきて，開始点・集結点を出力する。↑と被るからほぼdeprecated
def print_domain_gpdpos(alnfasta_dict):
    fa2_dict = alnfasta_dict

    for k,v in fa2_dict.items(): 
        sys.stderr.write("sequence length with gap: {}".format(len(v.seq)))
        break

    for domn in domain_name_dict.keys():
        dom_start,dom_end = 100000000000,0

        #先に最長となるようなDomainの開始・終了点を決めて，次のfor文で配列とってきて書き込んでいく
        for k_acc in fa2_dict.keys():
            try:
                dom_start_tmp = seqpos_to_gapped_seqpos(domain_dict[k_acc][domn]["start"],fa2_dict[k_acc].seq)
                dom_end_tmp   = seqpos_to_gapped_seqpos(domain_dict[k_acc][domn]["end"],fa2_dict[k_acc].seq)

                dom_start = dom_start_tmp if dom_start_tmp < dom_start else dom_start
                dom_end   = dom_end_tmp   if dom_end_tmp   > dom_end   else dom_end
            except KeyError as err: #主にそのドメインがヒットしなかった配列
                pass
        sys.stderr.write("{}\t{}\t{}".format(domn,dom_start,dom_end))
        
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

def make_domain_fasta(alnfasta_dict): #配列長，ドメインごとに揃えたりしない. アライメント後のfastaにする意味は特にない。
#strするときにpos-1で指定することを忘れずに！
    for domn in domain_name_dict.keys():
        for k_acc in alnfasta_dict.keys():
            try:
                dom_start = seqpos_to_gapped_seqpos(domain_dict[k_acc][domn]["start"],alnfasta_dict[k_acc].seq)
                dom_end   = seqpos_to_gapped_seqpos(domain_dict[k_acc][domn]["end"],alnfasta_dict[k_acc].seq)
            except KeyError as err: #主にそのドメインがヒットしなかった配列
                pass

            fasta_header = "{}-{}".format(k_acc, domn)
            sec_seq   = str(alnfasta_dict[k_acc].seq[dom_start-1:dom_end-1])
            sec_seq_nogap = sec_seq.replace('-','')
            outfa_fh.write(">{}\n{}\n".format(fasta_header,sec_seq_nogap))


##########
# main() #
##########
if __name__ == '__main__':

    #2016.3.8 HXB2
    outfa_fh = sys.stdout
    seq_sec = {'Integrase': [1437, 1744], 'p51_RT': [862, 1316], 'Protease': [762, 861], 'p15': [1317, 1436]}  #Integraseは 1728→1744へ。（最後尾まで含めるため）
    fa_d = read_fasta_todict("../data/removed-HXB2_mafft-linsi_with-HXB2_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta")
    make_seqsec_fasta(seq_sec, fa_d)
    quit()


    outfa_fh = sys.stdout
    read_hmmscan_domtblout_todict(hmmscan_file)
    fa_d = read_fasta_todict(in_alnfasta)
#    print_domain_gpdpos(fa_d)
    make_domain_fasta(fa_d)
    quit()
    
    seq_sec = {'prot_sec':[328,487], 'RT_sec':[488,873], 'RNase_sec':[874,1005], 'int_sec':[1006,1306]}
    seq_dom = get_domain_gpdpos_widest(fa_d)
    print(seq_dom)
    make_seqsec_fasta(seq_sec, fa_d)
    #Multiprocessing
    #pool = mp.Pool(proc)
    #cb = pool.map(make_fasta_by_domain,target_files)
    
    sys.stderr.write("{} [s]".format(time.time() - stime))
