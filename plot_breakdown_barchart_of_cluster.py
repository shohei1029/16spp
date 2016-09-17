
import time
stime = time.time()

import os
import sys
#import subprocess
#import multiprocessing as mp
import csv
from collections import defaultdict,OrderedDict
#from Bio import SeqIO
#import matplotlib.pyplot as plt

import pandas as pd

# Python -> R 連携用
import pyper as pr
r = pr.R(use_pandas = "True")
r("library(ggplot2)")
r("library(cowplot)")

argvs = sys.argv


###originalとの変更点等
#データセットの配列数で割り，割合をplotすることにする。

###その後の予定
#クラスター数で割って，データセットと同じ割合だとした場合の期待値も表示する。
#ex. クラスター数4だったら1/4とか
###

# 2016.8.4
# 15fallのEx_makecowplotR-Barchart_bycluster.py をベースに改造
# より新しいEx_makeR-Barchart_bycluster.pyをベースにしてないのは，そっちはデータセットの総数で割る作業を行っており，その部分がハードコーディングなため。

#############
# I/O files #
#############
in_file = argvs[1] #scpsの結果ファイルをinput?
core_withext = in_file.split('/')[-1]
core = os.path.splitext(core_withext)[0]
#original_fasta = "./fasta/HIV-1-gM-noRs_pol-aa_v3.fasta"
attrs_file = "../analysis/attrs_HIV,SIV_regionyear.txt" 
output_dir = "../results/cluster_breakdown/{}".format(core)

#make output dir
os.makedirs(output_dir,exist_ok=True)

print("Dealing with",in_file)
###########################
# Get data from scps file #
###########################
data_dict = defaultdict(lambda:defaultdict())
# 'acc': {'Region': 'Africa', 'Subtype': 'VER', 'ID': 'DJ048201|SIV|subtype:VER', 'cluster': '3', 'Organism_and_Subtype': 'SIV_VER', 'Region_and_Year': 'Africa_unknown', 'SamplingYear': '', 'Organism': 'SIV', 'Country': 'KENYA'}

cluster_set = set()
with open(in_file,'r') as in_fh:
    for line in in_fh:
        if "%" in line:
            continue
        if len(line) < 2: #空白行っぽいの飛ばす
            continue
        line = line.rstrip()
        header,cluster = line.split()

        data_dict[header]["Cluster"] = cluster
        cluster_set.add(cluster)


############################
# Get data from attrs file #
############################
### 解析に使われてない配列のAttributeデータも保存される。
with open(attrs_file,'r') as in_csv:
    for row in csv.DictReader(in_csv,delimiter = '\t'):
        data_dict[ row["ID"] ].update(row)
#        print(data_dict[row["ID"]])


##############
# Sum up ... #
##############
total_dict = defaultdict(lambda: defaultdict( lambda:defaultdict(int) ))
#{ 1:{subtype: {A:2,B:5,} , region: {america:500,} , .. },2:{hoge}
for c in cluster_set:
    for k,v in data_dict.items():
        try:
#            print(v["Cluster"])
            if v["Cluster"] == c:
                # logger
                #total_dict[c]["Organism"][ v["Organism"] ] += 1
                total_dict[c]["Subtype"][ v["Subtype"] ] += 1
                total_dict[c]["Country"][ v["Country"] ] += 1
                total_dict[c]["Region"][ v["Region"] ] += 1
                total_dict[c]["SamplingYear"][ v["SamplingYear"] ] += 1
                #total_dict[c]["Organism_and_Subtype"][ v["Organism_and_Subtype"] ] += 1
        except KeyError:
            pass

#print(total_dict["1"])
    
    
################################ 
# Calculate and draw pie chart #    
################################ 
def dict_to_list(dic):
    key_list = []
    value_list = []
    sorted_dic = OrderedDict(sorted(dic.items(), key=lambda t:t[0])) #Pie chartとBar chart ではソートに使うキーが異なる
#    print(sorted_dic)
    for k,v in sorted_dic.items():
        if k == '': # labelが''のものをとばす
            continue
        key_list.append(k)
        value_list.append(v)
    return key_list,value_list
        

print("Number of clusters:",len(total_dict.keys()))
for c,vdicts in total_dict.items():
    #make output dir
    output_dir2 = "{}/{}".format(output_dir,c)
    os.makedirs(output_dir2,exist_ok=True)

    for kkind,vdict in vdicts.items(): #result -> Subtype defaultdict(<class 'int'>, {'C': 380})
        #kkind : Subtype, SamplingYear, Region, Country
        label,data = dict_to_list(vdict)
        df = pd.DataFrame({'label':label, 'value':data})
        #print(df)

        # R - ggplot2
        r.assign('rdata',df)

#        save_name = "{}/{}.png".format(output_dir2,kkind)
        save_name = "{}/3plots.png".format(output_dir2)
        #p_title = "Cluster No.{}, {}".format(c,kkind)
        p_title = ""
        r.assign('save_name',save_name)
        r.assign('p_title',p_title)

#        print(r("rdata"))
        if kkind == "Subtype":
            r("plot.st <- ggplot(aes(x=label,y=value),data=rdata) + geom_bar(width=0.8,stat = 'identity') + labs(x='{xlab}', y='{ylab}', title=p_title)".format(xlab=kkind,ylab="Number of sequences"))
            r("plot.st <- plot.st + theme(axis.text.x = element_text(size=18,angle=0), \
                                          axis.text.y = element_text(size=14), \
                                          axis.title = element_text(size=18) )") #vjust=1で上揃え,hjustで文字列を揃えるのか！

        elif kkind == "SamplingYear":
            r("plot.sy <- ggplot(aes(x=label,y=value),data=rdata) + geom_bar(width=0.8,stat = 'identity') + labs(x='{xlab}', y='{ylab}', title=p_title)".format(xlab=kkind,ylab="Number of sequences"))
            r("plot.sy <- plot.sy + theme(axis.text.x = element_text(size=14,angle=70,hjust=1), \
                                          axis.text.y = element_text(size=14), \
                                          axis.title = element_text(size=18) )") #vjust=1で上揃え,hjustで文字列を揃えるのか！

        elif kkind == "Region":
            r("plot.rg <- ggplot(aes(x=label,y=value),data=rdata) + geom_bar(width=0.8,stat = 'identity') + labs(x='{xlab}', y='{ylab}', title=p_title)".format(xlab=kkind,ylab="Number of sequences"))
            r("plot.rg <- plot.rg + theme(axis.text.x = element_text(size=14,angle=70,hjust=1), \
                                          axis.text.y = element_text(size=14), \
                                          axis.title = element_text(size=18) )") #vjust=1で上揃え,hjustで文字列を揃えるのか！


#        r("save_plot(save_name,p,base_aspect_ratio = 1)")
#        r("ggsave(save_name,plot=p,width=7,height=7)")
    r("plots <- plot_grid(plot.st,plot.sy,plot.rg,labels = c('A','B','C'),ncol = 3,align = 'h')")
    r("save_plot(save_name,plots,base_aspect_ratio = 3)")
    print(r("save_name"),"正常に描画できたかはわからへんよ")
    



print(time.time()-stime,"[s]")
