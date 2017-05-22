#!/bin/env zsh
set -euo pipefail

# 2016.4.25 by Shohei Nagata
# 引数にname core(結果ファイル保存場所に使用)と，対象fastaを入れたらmkblastdbからblastやってくれるscript.

# 2016.8.18 memo
# v1(無印)との違いは，domains以外のregionにも対応したこと。引数で何のタイプか指定する。
# また，makeblastdb のオプションがエラー吐いたので，エラーをはかないように修正

# 2017.5.22 v3
# domainsとかregionとかPfamとかそういうハードコーディンの部分はなくした
# 引数は対象とするFASTAファイルのみに
# working ディレクトリ内にblastdbを作成。そこにPATHを通す方式へ

# sim作成スクリプトpythonのやつ，同じディレクトリ内におくように。暫定仕様。

#memo
#irisでやる設計
#stderrをログへ出力する形で実行したいね

#引数
# NAME CORE ex. HIV-1-gM-A,B,C-noRs_env-aa
# DOMAIN ex. RVT_1

EVALUE=1e-5
SCRIPT_BLAST_TO_SIM=../v5_blast7_tosims.py

FASTA_FILE=$1

WD_S=blastp_sim_$(date "+%Y%m%d")

mkdir -p ${WD_S}
cd ${WD_S}
WD=`pwd`

if [ $# -ne 1 ]; then
  echo "指定された引数は$#個です。" 1>&2
  echo "実行するには1個の引数が必要です。" 1>&2
  exit 1
fi

#########
# BLAST #
#########

#mkblastdb -> blastp @ iris
BLASTDB=${WD}/blastdb
export BLASTDB=${WD}/blastdb
mkdir -p ${BLASTDB}

if [[ ! -f ../${FASTA_FILE} ]]; then
    echo "can't find ${FASTA_FILE}, exitting.."
    exit 1
fi

#mkblastdb
echo "making BLAST Database.."
makeblastdb -in ../${FASTA_FILE} -out ${BLASTDB}/${FASTA_FILE} -dbtype prot -hash_index -parse_seqids -max_file_sz 2GB

#blast
BLAST_OUT_FILE_NAME=blastp_${EVALUE}.txt

echo "BLAST.."
blastp -query ../${FASTA_FILE} -db ${FASTA_FILE} -num_threads 16 -outfmt 7 -evalue ${EVALUE} -max_target_seqs 99999999 -out ${BLAST_OUT_FILE_NAME}


#######################
# create sif/sim file #
#######################

grep -v '#' ${BLAST_OUT_FILE_NAME} | LC_ALL=C sort -k 1,2 -u > mbs_${BLAST_OUT_FILE_NAME}

echo "creating sim files.."
cat mbs_${BLAST_OUT_FILE_NAME} | python ${SCRIPT_BLAST_TO_SIM} > sim_${BLAST_OUT_FILE_NAME}


# #%identityの閾値を変える時。未変更 #hoge.sim.txtじゃなく，hoge.txtになってしまう。。
# for pi in 0 {5..6}{0,5} {7..9}{0,2,4,5,6,8}
# do
#     cat ../analysis/blast/${DIR_NAME_CORE}/mbs_${BLAST_OUT_FILE_NAME} | python ./v5_blast7_tosims.py -p $pi > ../analysis/sim/${DIR_NAME_CORE}/pi${pi}_${BLAST_OUT_FILE_NAME}
#     echo $pi
# done
#
