
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
r = pr.R(use_pandas = True)
r("library(ggplot2)")
r("library(cowplot)")
r("library(data.table)")

argvs = sys.argv

#Options
no_xlabs_flag = True
no_title_flag = True



###今のところ，./Ex_makecowplotR-Barchart_bycluster.pyよりこっちがメイン
###originalとの変更点等
#データセットの配列数で割り，割合をplotすることにする。

#クラスター数で割って，データセットと同じ割合だとした場合の期待値も表示する。
#ex. クラスター数4だったら1/4とか

#scpsの結果ファイルを読み込ませる

#2015/11/24
#v2
#横軸を揃える。配列がなかったら0と表示。

#予定
#x軸なしオプションを付ける


#############
# I/O files #
#############
in_file = argvs[1] #scpsの結果ファイルをinput?
core_withext = in_file.split('/')[-1]
core = os.path.splitext(core_withext)[0]
#original_fasta = "./fasta/HIV-1-gM-noRs_pol-aa_v3.fasta"
#attrs_file = "../SCPS/attrs_HIV,SIV_regionyear.txt"
attrs_file = "../SCPS/attrs_domaincluster_std.txt" #for dommain attrs(IDにドメイン情報もついてるやつ)
if no_title_flag:
    output_dir = "../results/barchart_R_Ex_no-title/{}".format(core)
else:
    output_dir = "../results/barchart_R_Ex/{}".format(core)

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

#データセットの総数で割る
#ハードコーディング
dstotal_subtype = {"A1":142, "A2":2, "B":1528, "C":380}
dstotal_region  = {"Africa":413, "Asia":288, "Central_and_South_America":155, "Europe":185, "North_America":986, "Oceania":24}
dstotal_samplingyear ={ 
"1982":  1,
"1983": 10,
"1984":  4,
"1985":  8,
"1986": 14,
"1987":  2,
"1988":  1,
"1989": 38,
"1990": 15,
"1991": 12,
"1992": 14,
"1993": 12,
"1994":  8,
"1995": 50,
"1996": 74,
"1997": 21,
"1998": 47,
"1999": 45,
"2000":113,
"2001": 72,
"2002":103,
"2003":173,
"2004":232,
"2005":381,
"2006":229,
"2007":159,
"2008": 70,
"2009": 43,
"2010":  3,
"2011":  2}

dstotal_country = {
"RGENTINA":10,
"AUSTRALIA":24,
"BOLIVIA":1,
"BOTSWANA":38,
"BRAZIL":63,
"CAMEROON":1,
"CANADA":63,
"CHINA":46,
"COLOMBIA":2,
"CUBA":2,
"CYPRUS":97,
"DENMARK":17,
"DJIBOUTI":2,
"DOMINICANREPUBLIC":2,
"ETHIOPIA":1,
"FRANCE":9,
"GEORGIA":5,
"GERMANY":16,
"HAITI":5,
"HONGKONG":3,
"INDIA":16,
"ISRAEL":1,
"ITALY":1,
"JAMAICA":1,
"JAPAN":66,
"KAZAKHSTAN":6,
"KENYA":47,
"MALAWI":1,
"NETHERLANDS":12,
"PARAGUAY":3,
"PERU":58,
"RUSSIANFEDERATION":17,
"RWANDA":9,
"SENEGAL":4,
"SOUTHAFRICA":240,
"SOUTHKOREA":17,
"SPAIN":76,
"SWEDEN":3,
"SWITZERLAND":11,
"TANZANIA":22,
"THAILAND":22,
"TRINIDADANDTOBAGO":5,
"UGANDA":11,
"UKRAINE":13,
"UNITEDKINGDOM":10,
"UNITEDSTATES":923,
"URUGUAY":3,
"UZBEKISTAN":6,
"YEMEN":3,
"ZAMBIA":37}

#regions_set = {"Africa", "Asia", "Central_and_South_America", "Europe", "North_America", "Oceania"}
for cw in cluster_set:
    for kw,vw in total_dict[cw].items(): #ex. kw -> "Coutry"
        for ksp,vval in vw.items(): #ksp -> "JAPAN"
            try:
                if kw == "Subtype":
                    total_dict[cw]["Subtype"][ksp] = vval/dstotal_subtype[ksp]
                elif kw == "Country":
                    total_dict[cw]["Country"][ksp] = vval/dstotal_country[ksp]
                elif kw == "Region":
                    total_dict[cw]["Region"][ksp] = vval/dstotal_region[ksp]
                elif kw == "SamplingYear":
                    total_dict[cw]["SamplingYear"][ksp] = vval/dstotal_samplingyear[ksp]
            except KeyError as err:
                pass

##配列数0も0とするため
#for cw in cluster_set:
##    for kw in total_dict[cw].keys(): #ex. kw -> "Coutry"
#    for kw in list(total_dict[cw]): #ex. kw -> "Coutry"
##        for ksp,vval in vw.items():
#        if kw == "Region":
#            for r in regions_set:
#                if not total_dict[cw]["Region"][r]:
#                    total_dict[cw]["Region"][r] = 0.0
#                    print(total_dict[cw]["Region"][r])
#        elif kw == "SamplingYear":
#            for y in range(1982,2012):
#                if not total_dict[cw]["SamplingYear"][str(y)]:
#                    total_dict[cw]["SamplingYear"][str(y)] = 0.0
 
#期待値表示
expval = 1/len(cluster_set)
print("Expected value: ", expval)
 

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
#        print(kkind,vdict)
        label,data = dict_to_list(vdict)
        df = pd.DataFrame({'label':label, 'value':data})
        print(df)

        # R - ggplot2
        r.assign('rdata',df)
#        r("rdata = data.table(rdata)")

        save_name = "{}/{}.png".format(output_dir2,kkind)
        p_title = "Cluster No.{}, {}".format(c,kkind)
        r.assign('save_name',save_name)
        r.assign('p_title',p_title)

#        print(r("rdata"))
        r("p <- ggplot(aes(x=label,y=value),data=rdata) + geom_bar(width=0.8,stat = 'identity') + labs(x='{xlab}', y='{ylab}', title=p_title) + theme(axis.text.x = element_text(angle=70,)) + ylim(0,1) + geom_hline(aes(yintercept={yline}, colour='orange'),size=1)".format(xlab=kkind,ylab="Number of sequences",yline=expval))
        if kkind == "Subtype":
            r("p <- p + theme(axis.text.x = element_text(size=20,angle=0)) ") #vjust=1で上揃え,hjustで文字列を揃えるのか！
            if no_title_flag:
                r("p <- p + theme(title = element_blank())")
        elif kkind == "Region":
            r('p <- p + theme(axis.text.x = element_text(size=14,angle=70,hjust=1) ) + scale_x_discrete(limits=c("Africa", "Asia", "Central_and_South_America", "Europe", "North_America", "Oceania"), label=c("Central_and_South_America"="Central and\nSouthAmerica", "North_America"="North America"))') #vjust=1で上揃え,hjustで文字列を揃えるのか！
            #r('p <- p + theme(axis.text.x = element_text(size=14,angle=70,hjust=1) ) + xlim("Africa", "Asia", "Central_and_South_America", "Europe", "North_America", "Oceania") + scale_x_discrete(label=c("Central_and_South_America"="Central and\nSouthAmerica", "North_America"="North America"))') #vjust=1で上揃え,hjustで文字列を揃えるのか！
            #r('p <- p + theme(axis.text.x = element_text(size=14,angle=70,hjust=1) ) + xlim("Africa", "Asia", "Central_and_South_America", "Europe", "North_America", "Oceania")') #vjust=1で上揃え,hjustで文字列を揃えるのか！
            if no_title_flag:
                r("p <- p + theme(title = element_blank())")
            if no_xlabs_flag:
                r("p <- p + theme(axis.text.x = element_blank())")
        elif kkind == "SamplingYear":
            r("p <- p + theme(axis.text.x = element_text(size=14,angle=70,hjust=1) ) + xlim(paste(1982:2011,sep=','))") #vjust=1で上揃え,hjustで文字列を揃えるのか！
            if no_title_flag:
                r("p <- p + theme(title = element_blank())")
            if no_xlabs_flag:
                r("p <- p + theme(axis.text.x = element_blank())")
        else:
            r("p <- p + theme(axis.text.x = element_text(size=14,angle=70,hjust=1) ) ") #vjust=1で上揃え,hjustで文字列を揃えるのか！
            if no_title_flag:
                r("p <- p + theme(title = element_blank(), text = element_blank())")
#        r("save_plot(save_name,p,base_aspect_ratio = 1)")
        r("ggsave(save_name,plot=p,width=7,height=7)")
        print(r("save_name"))

    



print(time.time()-stime,"[s]")
