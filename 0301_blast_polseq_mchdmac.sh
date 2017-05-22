#!/opt/local/bin/zsh
set -euo pipefail

#mkblastdb -> blastp @ CNS (ex. mchdmac)

NAME_CORE=pol-4secs_A,B,C_HIV-1-gM-noRs_pol-aa_v3
BLASTDB=/home/t14650sn/bio/blastdb/

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

jobdone mkblastdb ${NAME_CORE} mchdmac

#blast
cd $WD

OUTPUT_FILE=../analysis/blast/${NAME_CORE}.txt

mkdir -p ../analysis/blast

blastp -query ${INPUT_FILE} -db ${NAME_CORE}.fasta -num_threads 16 -outfmt 7 -evalue 1e-5 -max_target_seqs 99999999 -out ${OUTPUT_FILE}

#end
jobdone blastp ${NAME_CORE} mchdmac
