#!/bin/env zsh
set -euo pipefail

# 2016.4.25 by Shohei Nagata
# 引数にname core(結果ファイル保存場所に使用)と，対象fastaを入れたらmkblastdbからblastやってくれるscript.

#memo
#irisでやる設計
#stderrをログへ出力する形で実行したいね

#引数
# NAME CORE ex. HIV-1-gM-A,B,C-noRs_env-aa
# DOMAIN ex. RVT_1

#ボツ
## FASTA FILE ex. ../data/domains_${NAME_CORE}/mafft_linsi_${DOMAIN}_domains-HXB2-Pfam_${NAME_CORE}.fasta # 変数は展開済みで入れる.

NAME_CORE=$1
DOMAIN=$2
#FASTA_FILE=$3
DIR_NAME_CORE="domains_${NAME_CORE}"
DETAILED_NAME_CORE="${DOMAIN}_domains-HXB2-Pfam_${NAME_CORE}"
FASTA_FILE="../data/${DIR_NAME_CORE}/${DETAILED_NAME_CORE}.fasta"

WD=`pwd`

if [ $# -ne 2 ]; then
  echo "指定された引数は$#個です。" 1>&2
  echo "実行するには3個の引数が必要です。" 1>&2
  exit 1
fi

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
makeblastdb -in ${DETAILED_NAME_CORE}.fasta -dbtype prot -hash_index -parse_seqids -max_file_sz 9000GB

#notify
jobdone mkblastdb ${FASTA_FILE} iris

#blast
cd $WD

BLAST_OUT_IRIS=/work/shohei1029/16spp/analysis/blast 
BLAST_OUT_FILE_NAME=blastp_1e-5_${DETAILED_NAME_CORE}.txt
OUTPUT_FILE=${BLAST_OUT_IRIS}/${BLAST_OUT_FILE_NAME}

mkdir -p $BLAST_OUT_IRIS #for iris
mkdir -p ../analysis/blast/${DIR_NAME_CORE}

cp $FASTA_FILE $BLAST_OUT_IRIS
cd $BLAST_OUT_IRIS
echo "BLAST.."
blastp -query ${DETAILED_NAME_CORE}.fasta -db ${DETAILED_NAME_CORE}.fasta -num_threads 72 -outfmt 7 -evalue 1e-5 -max_target_seqs 99999999 -out ${OUTPUT_FILE}

cp ${OUTPUT_FILE} ${WD}/../analysis/blast/${DIR_NAME_CORE}
cd $WD

#notify
jobdone blastp ${DETAILED_NAME_CORE} iris


#######################
# create sif/sim file #
#######################

mkdir -p ../analysis/sim/${DIR_NAME_CORE}
grep -v '#' ../analysis/blast/${DIR_NAME_CORE}/${BLAST_OUT_FILE_NAME} | LC_ALL=C sort -k 1,2 -u > ../analysis/blast/${DIR_NAME_CORE}/mbs_${BLAST_OUT_FILE_NAME}

echo "creating sim files.."
#hoge.sim.txtじゃなく，hoge.txtになってしまう。。
for pi in 0 {5..6}{0,5} {7..9}{0,2,4,5,6,8}
do
    cat ../analysis/blast/${DIR_NAME_CORE}/mbs_${BLAST_OUT_FILE_NAME} | python ./v5_blast7_tosims.py -p $pi > ../analysis/sim/${DIR_NAME_CORE}/pi${pi}_${BLAST_OUT_FILE_NAME}
    echo $pi
done

jobdone create sif/sim file ${DETAILED_NAME_CORE} iris
