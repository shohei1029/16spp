#!/usr/bin/env python3

import re
import time
import sys
#import numpy as np
from collections import defaultdict
import argparse
import logging

from progressbar import ProgressBar

logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
logger.setLevel(logging.INFO)
logger.addHandler(handler)

start_time = time.time()

# v3との変更点(change log)
# outpupt filename の変更
#  -> edge絞りの情報削除
# default で bothfmt
# --nomultihit も default で True (bit scoreの最大値のみ記録されているファイルを読み込ませる)

# * 類似度スコアでのcutoffを"未満"から"以下"へ変更 ***

# v4からの変更点(change log)
# scps用 sim 形式のみの出力へ (cytoscapeもその形式で読み込むし) <ID1 ID2 sim>
# ファイル出力時の拡張子変更 .sim から .sim.txtへ
# loggingの導入

#予定?
# defaultで，sys.stdinから読み込んでsys.stdoutへ出力 (unix!)


###################
# Argument Parser #
###################
parser = argparse.ArgumentParser()
parser.add_argument("-m","--nomultihit",type=bool,default=True,choices=[True,False],help="If there are NO multiple hits in Blast output file -> True .(分割ヒット対策をしなくていいんだったらTrue, 分割ヒットしててbit scoreの最大値を選んでとってくる必要があるんだったらFalse (Defauls)")
parser.add_argument("-i","--input_file",type=str,help="Input Blast output file (tabular).")
parser.add_argument("-w","--output_dir",type=str,help="Output directory")
parser.add_argument("-k","--edge_num",type=int,default=0,help="Ideal number of edges. Conserve alledges -> 0（エッジ数絞らなくていいなら0を指定")
parser.add_argument("-c","--cutoff_edge_weight",type=float,default=0.0,help="Cutoff by edge weight (i.e. sim(x,y)).")
args = parser.parse_args()

edge_num = args.edge_num
edge_weight = args.cutoff_edge_weight


#############
# I/O files #
#############
logger.info("1st phase ---> Setting I/O files")
in_file = args.input_file
out_dir = args.output_dir

if args.input_file:
    in_fh = open(in_file,'r'); logger.info('1st phase ---> Opening...' + in_file) #BLAST output file (tabular format)
    out_fh = sys.stdout
else:
    in_fh = sys.stdin
    
if args.output_dir:
    if edge_num == 0:
        out_file = '{}/t{}_{}.sim.txt'.format(out_dir,edge_weight,in_file.split('/')[-1])
    else:
        logger.info("edge SHIBORIKOMI mode")
        out_file = '{}/edges{}_t{}_{}.sim.txt'.format(out_dir,edge_num,edge_weight,in_file.split('/')[-1])
    out_fh = open(out_file,'w')
else:
    out_fh = sys.stdout
    out_file = "sys.stdout"


###########################
# Parse Blast output file #
###########################
#0->query,1->subject,11->bit scores
logger.info("2nd phase ---> Parse Blast output file...")
blast = defaultdict(dict)
tmp_scores = defaultdict(dict)
entry_set = set()
p = re.compile('#')

if args.nomultihit == True:
    logger.info("NO Multi hit Blast ! :D")
    for line in in_fh:
        if p.match(line):
            continue
        line = line.rstrip()
    
        query   = line.split('\t')[0]
        subject = line.split('\t')[1]
        score   = float(line.split('\t')[11])

        blast[query][subject] = score
    
        entry_set.add(query)
        entry_set.add(subject)
    
else:
    logger.info("Multi hit Blast mode !!")
    for line in in_fh:
        if p.match(line):
            continue
        line = line.rstrip()
    
        query   = line.split('\t')[0]
        subject = line.split('\t')[1]
        score   = float(line.split('\t')[11])
    
        #分割ヒット対策をしたい。最も大きなbit scoreの値のみを入れたい
        try:
            tmp_scores[query][subject].append(score) 
        except KeyError as err:
            tmp_scores[query][subject] = []
            tmp_scores[query][subject].append(score) 
    
        blast[query][subject] = float(max(tmp_scores[query][subject]))
    
        entry_set.add(query)
        entry_set.add(subject)
    

#########################################################
# Calculate sequence similarity scores (Matsui metrics) #
#########################################################
logger.info("3rd phase ---> Calculate sequence similarity scores")
entry_list = list(entry_set)
entry_list = sorted(entry_list) #not necessary?
mat = defaultdict(dict)

pbar = ProgressBar(max_value=len(entry_list),redirect_stderr=True)
pbar.start()

#entry_list -> ID str
#blast[str][str]
#mat[str][str]
othererr = []
pbar_count = 0
for i in entry_list:
    for j in entry_list:
        ij=0.0; ji=0.0; ii=0.0; jj=0.0;

        try:
            ij = blast[i][j]
        except KeyError:
            ij = 0
        try:
            ji = blast[j][i]
        except KeyError:
            ji = 0
        try:
            ii = blast[i][i]
        except KeyError:
            ii = 0
        try:
            jj = blast[j][j]
        except KeyError:
            jj = 0

        dist  = ij if ij >= ji else ji
        self  = ii if ii >= jj else jj #self == ii としたら初期スクリプトと同じになるかなぁ→でもこっちのほうが論文にも載ってるし正しい

        score = 0.0
        try:
            score = dist/self
        except ZeroDivisionError as err:
            score = dist/1
            logger.warn("ZeroDivisionError!" + err)
            othererr.append(err)

        mat[i][j] = 0 if i == j else score #松井さんのPerlスクリプトでは0にするようになっていたが，やめる。→どっちでもいいん？
#        mat[i][j] = score #どっちでもいいみたい？でも結果なぁんか違うような気もする。その程度とも言える。
    pbar_count += 1
    pbar.update(pbar_count)
pbar.finish()


###############################
# Select top 5 edges per node #
###############################
logger.info("4th phase ---> Select top {} edges per node".format(edge_num))
best = defaultdict(dict)
#best[str][str]

if edge_num == 0:
    logger.info("Conserve all edges.")
    best = mat
else: #1つのノードと接続しているエッジ数を限定する場合
    #行列の場所をIDでやる。IDの場所を別のリストと対応づけてその場所の番号（数字）で比較するようにする。
    for i,ic in enumerate(entry_list):
        sorted_list = sorted(mat[ic], reverse = True) #降順（大きい順）
        #print(sorted_list)
        for j,jc in enumerate(sorted_list[0:edge_num]):
    
            l,r = "",""
            if i>= j: #対称行列の半分しか出力しないため，出力する側へ移している
                l,r = ic,jc
            else:
                l,r = jc,ic
                
            best[l][r] = mat[ic][jc]

            #暫定 (でも他に思いつかない)
            if i >= j and best[ic][jc] > edge_weight:
                out_fh.write("{} {} {}\n".format(i,best[i][j],j)) #sif
    logger.info("暫定システム")
    out_fh.close()
    quit()


#######################################################
# Output network as cytoscape or/and SCPS format file #
#######################################################
logger.info("5th phase ---> Output network as Cytoscape or/and SCPS format file")

#best[str][str]
for i in best.keys():
    for j in best[i].keys():
        if i >= j and best[i][j] > edge_weight:
            out_fh.write("{} {} {}\n".format(i,best[i][j],j)) #sif


##########
# finish #
##########
logger.info("End Phase ---><")
out_fh.close()
in_fh.close()

logger.info("Output file -> " + out_file)
logger.info("Run time: "+ str((time.time() - start_time)/60) + '[m]')


############
# for test #
############
#test_fh = open("./test_matrix_self1_hiv2norecombinsnts.csv","w")
#for i in best.keys():
#    logger.info(i)
#    for j in best[i].keys():
#        test_fh.write("{},".format(best[i][j]))
#    test_fh.write("\n")
#
#test_fh.close()
