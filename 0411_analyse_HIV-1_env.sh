#!/bin/env zsh
set -euo pipefail

# 2016.4.11 by Shohei Nagata
# env領域のアミノ酸配列抽出，総当りblast，sim作成，mafftによるアライメントまで
# v1.0.1
#  生成されるsimファイル名のバグを修正

#memo
#irisでやる設計
#とりあえず，サブタイプをA, B, Cに絞らずやってみる
#stderrをログへ出力する形で実行したいね

NAME_CORE=HIV-1-gM-noRs_env-aa
FASTA_FILE="../data/${NAME_CORE}.fasta"
WD=`pwd`

########################
# make data sets #
########################
#-t:target gene/protein name , -m s: specify subtype (remove recombinants) mode , -s: specify subtype (sep=,)| -n: remove gene name in fasta header, criteria default: 90%, refseqlength: isolate HXB2
python getseq_ingbk_withbginfo_tofasta.py -t env -m s -s A1,A2,B,C,D,F1,F2,G,H,J,K | python filter_fasta_by_seq-length.py -n True -r 856 > $FASTA_FILE

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
makeblastdb -in ${NAME_CORE}.fasta -dbtype prot -hash_index -parse_seqids -max_file_sz 9000GB

#notify
jobdone mkblastdb ${NAME_CORE} iris

#blast
cd $WD

BLAST_OUT_IRIS=/work/shohei1029/16spp/analysis/blast 
BLAST_OUT_FILE_NAME=blastp_1e-5_${NAME_CORE}.txt
OUTPUT_FILE=${BLAST_OUT_IRIS}/${BLAST_OUT_FILE_NAME}

mkdir -p $BLAST_OUT_IRIS #for iris
mkdir -p ../analysis/blast

blastp -query ${FASTA_FILE} -db ${NAME_CORE}.fasta -num_threads 72 -outfmt 7 -evalue 1e-5 -max_target_seqs 99999999 -out ${OUTPUT_FILE}

cp ${OUTPUT_FILE} ../analysis/blast

#notify
jobdone blastp ${NAME_CORE} iris


#######################
# create sif/sim file #
#######################

mkdir -p ../analysis/sim/${NAME_CORE}
grep -v '#' ../analysis/blast/${BLAST_OUT_FILE_NAME} | LC_ALL=C sort -k 1,2 -u > ../analysis/blast/mbs_${BLAST_OUT_FILE_NAME}

#hoge.sim.txtじゃなく，hoge.txtになってしまう。。
for pi in 0 60 70 75 80 85 90 95
do
    cat ../analysis/blast/mbs_${BLAST_OUT_FILE_NAME} | python ./v5_blast7_tosims.py -p $pi > ../analysis/sim/${NAME_CORE}/pi${pi}_${BLAST_OUT_FILE_NAME}
done

jobdone create sif/sim file ${NAME_CORE} iris

#############
# alignment #
#############
mafft-linsi --thread 72 $FASTA_FILE > ../data/mafft-linsi_${NAME_CORE}.fasta

jobdone mafft `mafft --version` ${NAME_CORE} iris


# transfer files #
rsync -av ../analysis/sim/${NAME_CORE} shohei@133.27.17.109:/home/shohei/bio/16spp/analysis/sim
