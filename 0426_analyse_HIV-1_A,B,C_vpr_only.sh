#!/bin/env zsh
set -euo pipefail

# 2016.4.11 by Shohei Nagata
# env領域のアミノ酸配列抽出，総当りblast，sim作成，mafftによるアライメントまで
# v1.0.1
#  生成されるsimファイル名のバグを修正
# 4.26
# vif 用へ。ドメイン等なしで。

#memo
#irisでやる設計
#とりあえず，サブタイプをA, B, Cに絞らずやってみる
#stderrをログへ出力する形で実行したいね

NAME_CORE=HIV-1-gM-A,B,C-noRs_vpr-aa
FASTA_FILE_NAME=${NAME_CORE}.fasta
FASTA_FILE="../data/${FASTA_FILE_NAME}" 
PROT_LENGTH=96
TARGET_NAME=vpr
WD=`pwd`

########################
# make data sets #
########################
#-t:target gene/protein name , -m s: specify subtype (remove recombinants) mode , -s: specify subtype (sep=,)| -n: remove gene name in fasta header, criteria default: 90%, refseqlength: isolate HXB2
python getseq_ingbk_withbginfo_tofasta.py -t ${TARGET_NAME} -m s -s A1,A2,B,C | python filter_fasta_by_seq-length.py -n True -r ${PROT_LENGTH} > $FASTA_FILE
echo 'number of sequences: '
cat $FASTA_FILE | grep '>' | wc -l

#notify
jobdone make data sets ${NAME_CORE} iris


#########
# BLAST #
#########

#mkblastdb -> blastp @ iris
BLASTDB=${HOME}/blastdb

if [[ ! -f ${FASTA_FILE} ]]; then
    echo "can't find ${FASTA_FILE}, exitting.."
    exit 1
fi

#mkblastdb
cp ${FASTA_FILE} ${BLASTDB}
cd $BLASTDB
makeblastdb -in ${FASTA_FILE_NAME} -dbtype prot -hash_index -parse_seqids -max_file_sz 9000GB

#notify
jobdone mkblastdb ${NAME_CORE} iris

#blast
cd $WD

BLAST_OUT_IRIS=/work/shohei1029/16spp/analysis/blast 
BLAST_OUT_FILE_NAME=blastp_1e-5_${NAME_CORE}.txt
OUTPUT_FILE=${BLAST_OUT_IRIS}/${BLAST_OUT_FILE_NAME}

mkdir -p $BLAST_OUT_IRIS #for iris
mkdir -p ../analysis/blast

cp $FASTA_FILE $BLAST_OUT_IRIS
cd $BLAST_OUT_IRIS
echo "BLAST.."
blastp -query ${FASTA_FILE_NAME} -db ${FASTA_FILE_NAME} -num_threads 72 -outfmt 7 -evalue 1e-5 -max_target_seqs 99999999 -out ${OUTPUT_FILE}

cp ${OUTPUT_FILE} ${WD}/../analysis/blast
cd $WD

#notify
jobdone blastp ${NAME_CORE} iris


#######################
# create sif/sim file #
#######################

mkdir -p ../analysis/sim/${NAME_CORE}
grep -v '#' ../analysis/blast/${BLAST_OUT_FILE_NAME} | LC_ALL=C sort -k 1,2 -u > ../analysis/blast/mbs_${BLAST_OUT_FILE_NAME}

echo "creating sim files.."
#hoge.sim.txtじゃなく，hoge.txtになってしまう。。
for pi in 0 {5..6}{0,5} {7..9}{0,2,4,5,6,8}
do
    cat ../analysis/blast/mbs_${BLAST_OUT_FILE_NAME} | python ./v5_blast7_tosims.py -p $pi > ../analysis/sim/${NAME_CORE}/pi${pi}_${BLAST_OUT_FILE_NAME}
done

jobdone create sif/sim file ${NAME_CORE} iris

#############
# alignment #
#############
MAFFT_WD=/work/shohei1029/16spp/analysis/mafft
mkdir -p ${MAFFT_WD}
cp $FASTA_FILE $MAFFT_WD
cd $MAFFT_WD

mafft-linsi --thread 72 $FASTA_FILE_NAME > ${WD}/../data/mafft-linsi_${NAME_CORE}.fasta

cd $WD

jobdone mafft `mafft --version` ${NAME_CORE} iris


# transfer files #
rsync -av ../analysis/sim/${NAME_CORE} shohei@133.27.17.109:/home/shohei/bio/16spp/analysis/sim
