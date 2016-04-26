
import sys
import argparse
import os
import logging #log
from collections import defaultdict,Counter
from Bio import SeqIO
argv = sys.argv


logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setLevel(logging.WARN)
logger.setLevel(logging.WARN)
logger.addHandler(handler)

# なぜかSubtypes指定で，HIV-1とHIV-2の両方から配列をもってこれない！

# 2016.4.9
# copied from 14fall
# modified
# 出力fastaのheaderにgene nameがついているのはそのままにする。配列長filteringするスクリプトを作って，そこで遺伝子名除去もやる


###################
# Argument Parser #
###################
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--organism", default="HIV-1", type=str, help="HIV-1,HIV-2,SIV", choices=["HIV-1","HIV-2","SIV"])
#parser.add_argument("-t", "--target", type=str, help="What kind of data do you want to  get?", choices=["Patient_code","Accession","Subtype","Country","SamplingYear","Seq_Length"])
parser.add_argument("-t", "--target_gene", default="env", type=str, help="target gene/protein name")
parser.add_argument("-s", "--target_subtypes", type=str, help="target groups/subtypes. sep=, . ex. A1,A2,B,C ..")
parser.add_argument("-m", "--mode", default="n", choices=["c","s","r","n"], help="c:get seqs which have Country information, s:specify Subtypes")
args = parser.parse_args()
###################

fh = open("../data/data_sets/HIV,SIV_info.txt",'r')
output_fasta = "./test/test_HIV-1.fasta"

target_gene_aa = args.target_gene
if args.target_subtypes:
    target_subtypes = set(args.target_subtypes.split(","))
else:
    target_subtypes = set(["A1","A2","B","C","D","F1","F2","G","H","J","K","N","O","P","U","A"]) #HIV-1,HIV-2 except recombinants
#random_seq_num = 60


#############
# Functions #
#############
database_dic = {}
gb_acc = set() #条件に一致した配列のAccessionを記録。genbankファイルよりもってきてheader情報と合わせてFASTAへ
num_assersionerr = 0
def opengb_recursively(acc_list=gb_acc):
#    output_handle = open(output_fasta,'w') 
    output_handle = sys.stdout
    for root, dirs, files in os.walk('../data/data_sets/'):
        for file_ in files:
            if os.path.splitext(file_)[1] == '.genbank':
                logger.info(file_ )
                with open(os.path.join(root, file_)) as input_handle:
                    #get_origin_byacc_tofasta(input_handle,output_handle,acc_list) #切り替えできるようにしなきゃ...
                    get_geneaa_byacc_tofasta(input_handle,output_handle,acc_list)


def get_origin_byacc_tofasta(gbfile,output_fasta,acc_list):
    try:
        for seq_record in SeqIO.parse(gbfile, "genbank") :
            logger.info(seq_record.id)
#            if "genome" not in seq_record.description:
#                continue
            for a in acc_list:
                if a == seq_record.id.split('.')[0]:
                    logger.info("Dealing with GenBank record {}".format(seq_record.name))
                    output_fasta.write(">{acc}|{org}|subtype:{subtype}\n{seq}\n".format(
                        acc = seq_record.id.split('.')[0],
                        org = database_dic[seq_record.id.split('.')[0]]["Organism"],
                        subtype = database_dic[seq_record.id.split('.')[0]]["Subtype"],
                        seq = seq_record.seq))
    except AssertionError as err:
        num_assersionerr =+ 1
        logger.warn("AssertionError has been occured: " + str(err) )

no_target_aa = []
def get_geneaa_byacc_tofasta(gbfile,output_fasta,acc_list):
    try:
        for seq_record in SeqIO.parse(gbfile, "genbank") :
            logger.info(seq_record.id)
            for a in acc_list:
                if a == seq_record.id.split('.')[0]:
                    flag = False
                    logger.info("Dealing with GenBank record {}".format(seq_record.name))
                    for seq_feature in seq_record.features:
                        if seq_feature.type == 'CDS':
                            if "gene" in seq_feature.qualifiers:
                                if seq_feature.qualifiers["gene"][0] == target_gene_aa:
                                    logger.info("{} 'gene' found".format(target_gene_aa))
                                    flag = True
                                    try:
                                        output_fasta.write(">{acc}|{org}|subtype:{subtype}|{geneaa}\n{seq}\n".format(
                                        acc = seq_record.id.split('.')[0],
                                        org = database_dic[seq_record.id.split('.')[0]]["Organism"],
                                        subtype = database_dic[seq_record.id.split('.')[0]]["Subtype"],
                                        geneaa = seq_feature.qualifiers["gene"][0],
                                        seq = seq_feature.qualifiers["translation"][0]))
                                    except KeyError as err:
                                        logger.warn("No 'gene' translation!" + str(err) + " @" + seq_record.id)
                            elif "product" in seq_feature.qualifiers:
                                if seq_feature.qualifiers["product"][0] == "{} protein".format(target_gene_aa) or seq_feature.qualifiers["product"][0] == "{} polyprotein".format(target_gene_aa) :
                                    logger.info("{} 'product' found".format(target_gene_aa) + " @" + seq_record.id)
                                    flag = True
                                    try:
                                        output_fasta.write(">{acc}|{org}|subtype:{subtype}|{geneaa}\n{seq}\n".format(
                                        acc = seq_record.id.split('.')[0],
                                        org = database_dic[seq_record.id.split('.')[0]]["Organism"],
                                        subtype = database_dic[seq_record.id.split('.')[0]]["Subtype"],
                                        geneaa = seq_feature.qualifiers["product"][0],
                                        seq = seq_feature.qualifiers["translation"][0]))
                                    except KeyError as err:
                                        logger.warn("No 'product' translation!" + str(err) + " @" + seq_record.id)
                            else:
                                logger.info("No {} aa sequence in this CDS part.".format(target_gene_aa) + seq_record.id)
                    if flag == False: 
                        no_target_aa.append(a)
                        #なにかおかしい。ちゃんとpolがとれてきてるのにこのリストに入ってる奴もある。
    except AssertionError as err:
        num_assersionerr =+ 1
        logger.warn("AssertionError has been occured: " + str(err) )

def add_seqs_randomly_specify_subtypes(subtypes,rand=50):
    """subypes: dict{HIV-1:{subtype_name:{set}}, HIV-2:{subtype_name:{set}} }"""
    for dict_by_organims in subtypes.values():
        for v_set in dict_by_organims.values():
            leng = len(v_set)
            if leng > rand:
                leng = rand
            for i in range(leng):
                try:
                    gb_acc.add(v_set.pop())
                except KeyError:
                    break

#main()
########################
# Parce BackgroundInfo #
########################
for line in fh:
    line = line.rstrip()
    data = line.split('\t')
    #Accessionをkeyとしてそのvalueにいろんな情報(including Accession)が入ってる
    database_dic[data[3]] = {'Patient_code':data[1],'Accession':data[3],'Name':data[4],'Subtype':data[5],'Country':data[6],'SamplingYear':data[7],'Seq_Length':data[10],'Organism':data[11]}
########################


#条件に一致した配列のAccessionを記録。genbankファイルよりもってきてheader情報と合わせてFASTAへ
#Country 情報をもつもののみ
if args.mode == "c":
    logger.info("Get seqs which have Country information.")
#    if args.organism == "HIV-1" or args.organism == "HIV-2" or args.organism == "SIV":
    if args.organism: 
        gb_acc = [d["Accession"] for d in database_dic.values() if d['Organism'] == args.organism and d["Country"] != ""]
        opengb_recursively(gb_acc)
    else:
        gb_acc = [d["Accession"] for d in database_dic.values() if d["Country"] != ""]
        opengb_recursively(gb_acc)

#Subytpesを指定（うえのほう）
if args.mode == "s":
    logger.info("Specify Subtypes.")
    if args.organism: 
        logger.info("Specify Organims!")
        gb_acc = {d["Accession"] for d in database_dic.values() if d['Organism'] == args.organism and d["Subtype"] in target_subtypes}
        opengb_recursively(gb_acc)
    else:
#        gb_acc = [d["Accession"] for d in database_dic.values() if d["Subtype"] in target_subtypes and d['Organism'] != "FIV"]   #HIVのみを前提とする
        gb_acc = {d["Accession"] for d in database_dic.values() if d["Subtype"] in target_subtypes }   #HIVのみを前提とする
        opengb_recursively(gb_acc)

#多すぎるSbutypeからは50くらいランダムでもってくる
subtypes_source=defaultdict(lambda: defaultdict(set)) 
if args.mode == "r":
    logger.info("Specify Subtypes and choose {} seqs randumly.".format(random_seq_num))
    if args.organism: 
        logger.info("Specify Organims!")
        for d in database_dic.values():
            if d['Organism'] == args.organism and d["Subtype"] in target_subtypes:
                subtypes_source[d["Organism"]][d["Subtype"]].add(d["Accession"])

        add_seqs_randomly_specify_subtypes(subtypes_source,random_seq_num)
        opengb_recursively(gb_acc)

    else:
        for d in database_dic.values():
#            if d["Subtype"] in target_subtypes:
#                subtypes_source[d["Subtype"]].add(d["Accession"])
            if d['Organism'] != "SIV" and d["Subtype"] in target_subtypes: 
                subtypes_source[d["Organism"]][d["Subtype"]].add(d["Accession"])

        add_seqs_randomly_specify_subtypes(subtypes_source,random_seq_num)
        opengb_recursively(gb_acc)

#Country 情報の有無によらない
if args.mode == "n":
    logger.info("GN mode")
    gb_acc = {d["Accession"] for d in database_dic.values() if d['Organism'] != "FIV" } #HIVのみ
    opengb_recursively(gb_acc)


##########
# Finish #
##########
logger.info("{} times of AssertionError!".format(num_assersionerr))
logger.info("No {}: ".format(target_gene_aa) + str(no_target_aa) )
fh.close()


