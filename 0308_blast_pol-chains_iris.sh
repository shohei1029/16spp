#!/bin/env zsh
set -euo pipefail

#mkblastdb -> blastp @ iris

NAME_CORE=pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3
BLASTDB=${HOME}/blastdb

WD=`pwd`

INPUT_FILE="../data/${NAME_CORE}.fasta"

if [[ ! -f ${INPUT_FILE} ]]; then
    echo "can't find ${INPUT_FILE}, exitting.."
    exit 1
fi

#mkblastdb
cp ${INPUT_FILE} ${BLASTDB}
cd $BLASTDB
makeblastdb -in ${NAME_CORE}.fasta -dbtype prot -hash_index -parse_seqids -max_file_sz 9000GB

#notify
jobdone mkblastdb ${NAME_CORE} iris

#blast
cd $WD

BLAST_OUT_IRIS=/work/shohei1029/16spp/analysis/blast 
OUTPUT_FILE=${BLAST_OUT_IRIS}/blastp_1e-5_${NAME_CORE}.txt

mkdir -p $BLAST_OUT_IRIS #for iris
mkdir -p ../analysis/blast

blastp -query ${INPUT_FILE} -db ${NAME_CORE}.fasta -num_threads 72 -outfmt 7 -evalue 1e-5 -max_target_seqs 99999999 -out ${OUTPUT_FILE}

cp ${OUTPUT_FILE} ../analysis/blast


#end
jobdone blastp ${NAME_CORE} iris
