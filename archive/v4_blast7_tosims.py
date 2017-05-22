#!/usr/bin/env python3

import re
import time
#import sys
#import numpy as np
from collections import defaultdict
import argparse

start_time = time.time()

# v3との変更点(change log)
# outpupt filename の変更
#  -> edge絞りの情報削除
# default で bothfmt
# --nomultihit も default で True (bit scoreの最大値のみ記録されているファイルを読み込ませる)

# * 類似度スコアでのcutoffを"未満"から"以下"へ変更 ***


###################
# Argument Parser #
###################
parser = argparse.ArgumentParser()
parser.add_argument("-m","--nomultihit",type=bool,default=True,choices=[True,False],help="If there are NO multiple hits in Blast output file -> True .(分割ヒット対策をしなくていいんだったらTrue, 分割ヒットしててbit scoreの最大値を選んでとってくる必要があるんだったらFalse (Defauls)")
parser.add_argument("-i","--input_file",type=str,help="Input Blast output file (tabular).")
parser.add_argument("-k","--edge_num",type=int,default=0,help="Ideal number of edges. Conserve alledges -> 0（エッジ数絞らなくていいなら0を指定")
parser.add_argument("-s","--simfmt",type=bool,default=False,choices=[True,False],help="for SCPS format(ID1 ID2 sim)")
parser.add_argument("-b","--bothfmt",type=bool,default=True,choices=[True,False],help="Output both sim(SCPS) and sif(Cytoscape) format")
parser.add_argument("-c","--cutoff_edge_weight",type=float,default=0.0,help="Cutoff by edge weight (i.e. sim(x,y)).")
args = parser.parse_args()

edge_num = args.edge_num
edge_weight = args.cutoff_edge_weight


#############
# I/O files #
#############
in_file = args.input_file

if args.bothfmt == True:
    #out_file_sim = './sims/forscps_{}edges_sims_{}.sim'.format(edge_num,in_file.split('/')[-1])
    #out_file_sif = './sif/{}edges_sims_{}.sif'.format(edge_num,in_file.split('/')[-1])
    if edge_num == 0:
        out_file_sim = './sims/t{}_{}.sim'.format(edge_weight,in_file.split('/')[-1])
        out_file_sif = './sif/t{}_{}.sif'.format(edge_weight,in_file.split('/')[-1])
    else:
        print("edge SHIBORIKOMI mode")
        out_file_sim = './sims/edges{}_t{}_{}.sim'.format(edge_num,edge_weight,in_file.split('/')[-1])
        out_file_sif = './sif/edges{}_t{}_{}.sif'.format(edge_num,edge_weight,in_file.split('/')[-1])

elif args.simfmt == True:
    out_file = './sims/t{}_{}.sim'.format(edge_weight,in_file.split('/')[-1])
else:
    out_file = './sif/t{}_{}.sif'.format(edge_weight,in_file.split('/')[-1])

in_fh = open(in_file,'r'); print('1st phase ---> Opening...',in_file) #BLAST output file (tabular format)

if args.bothfmt == True:
    out_fh_sim = open(out_file_sim,'w')
    out_fh_sif = open(out_file_sif,'w')
else:
    out_fh = open(out_file,'w')

print(time.time() - start_time,"[s]")


###########################
# Parse Blast output file #
###########################
#0->query,1->subject,11->bit scores
print("2nd phase ---> Parse Blast output file...")
time_2nd = time.time()
blast = defaultdict(dict)
tmp_scores = defaultdict(dict)
entry_set = set()
p = re.compile('#')

if args.nomultihit == True:
    print("NO Multi hit Blast !!")
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
    print("Multi hit Blast mode !!")
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
    

print(time.time() - time_2nd,"[s]")


#########################################################
# Calculate sequence similarity scores (Matsui metrics) #
#########################################################
print("3rd phase ---> Calculate sequence similarity scores")
time_3rd = time.time()
entry_list = list(entry_set)
entry_list = sorted(entry_list) #not necessary?
mat = defaultdict(dict)

#entry_list -> ID str
#blast[str][str]
#mat[str][str]
othererr = []
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
            print("ZeroDivisionError!!!!!!!",err)
            othererr.append(err)

        mat[i][j] = 0 if i == j else score #松井さんのPerlスクリプトでは0にするようになっていたが，やめる。→どっちでもいいん？
#        mat[i][j] = score #どっちでもいいみたい？でも結果なぁんか違うような気もする。その程度とも言える。

print(time.time() - time_3rd,"[s]")


###############################
# Select top 5 edges per node #
###############################
print("4th phase ---> Select top {} edges per node".format(edge_num))
time_4th = time.time()
best = defaultdict(dict)
#best[str][str]

if edge_num == 0:
    print("Conserve all edges.")
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
                out_fh_sim.write("{} {} {}\n".format(l,r,best[l][r])) #sim for SCPS
                out_fh_sif.write("{} {} {}\n".format(l,best[l][r],r)) #sif for Cytoscape
    print("暫定システム")
    out_fh_sim.close()
    out_fh_sif.close()
    quit()


print(time.time() - time_4th,"[s]")


#######################################################
# Output network as cytoscape or/and SCPS format file #
#######################################################
print("5th phase ---> Output network as Cytoscape or/and SCPS format file")
time_5th = time.time()

#best[str][str]
if args.bothfmt == True:
    for i in best.keys():
        for j in best[i].keys():
            if i >= j and best[i][j] > edge_weight: #v4 未満から以下をカットするようへ変更
                out_fh_sim.write("{} {} {}\n".format(i,j,best[i][j])) #sim for SCPS
                out_fh_sif.write("{} {} {}\n".format(i,best[i][j],j)) #sif for Cytoscape
    out_fh_sim.close()
    out_fh_sif.close()
elif args.simfmt == True:
    for i in best.keys():
        for j in best[i].keys():
            if i >= j and best[i][j] > edge_weight:
                out_fh.write("{} {} {}\n".format(i,j,best[i][j])) #sim for SCPS
    out_fh.close()
else:
    for i in best.keys():
        for j in best[i].keys():
            if i >= j and best[i][j] > edge_weight:
                out_fh.write("{} {} {}\n".format(i,best[i][j],j)) #sif
    out_fh.close()

print(time.time() - time_5th,"[s]")


##########
# finish #
##########
print("End Phase ---><")
in_fh.close()

if args.bothfmt == True:
    print("Output file (sim for SCPS) ->",out_file_sim)
    print("Output file (sif for Cytoscape) ->",out_file_sif)
else:
    print("Output file ->",out_file)

print("Run time:",(time.time() - start_time)/60,'[m]')


############
# for test #
############
#test_fh = open("./test_matrix_self1_hiv2norecombinsnts.csv","w")
#for i in best.keys():
#    print(i)
#    for j in best[i].keys():
#        test_fh.write("{},".format(best[i][j]))
#    test_fh.write("\n")
#
#test_fh.close()
