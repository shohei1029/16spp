#!/bin/env zsh
set -euo pipefail

# 2016.8.9 by Shohei Nagata
# fastaファイルを読み込ませて，MAFFTでアライメントして，RAxMLで最尤系統樹構築

NAME_CORE="2seqs_eachcluster_scps-c5_blastp_1e-5_A,B,C_HIV-1-gM-noRs_pol-aa_v3_+OG-HIV-2"
FASTA_FILE="../data/${NAME_CORE}.fasta"
ALIGNED_FASTA="../data/mafft-linsi_${NAME_CORE}.fasta"
RAXML_OUTDIR="../analysis/RAxML/${NAME_CORE}"
WD=`pwd`

mkdir -p ${RAXML_OUTDIR}
RAXML_OUTDIR=`cd ${RAXML_OUTDIR} && pwd`
echo $RAXML_OUTDIR

if [[ ! -f ${FASTA_FILE} ]]; then
    echo "can't find ${FASTA_FILE}, exitting.."
    exit 1
fi


sed -i -e 's/:/_/g' ${FASTA_FILE}

mafft-linsi --thread 8 $FASTA_FILE > ${ALIGNED_FASTA}

raxmlHPC-PTHREADS-SSE3 -f a -p 12345 -x 12345 -m PROTGAMMAAUTO -N 1000 -o "D00835|HIV-2|subtype_A" -s ${ALIGNED_FASTA} -w ${RAXML_OUTDIR} -n ${NAME_CORE} -T 8 
