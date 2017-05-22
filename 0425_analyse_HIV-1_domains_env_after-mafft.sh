#!/bin/env zsh
set -euo pipefail

# 2016.4.24 by Shohei Nagata
# ドメインごとの領域を対象にsimファイルの生成まで。
# ドメインはisolate HXB2を基準に。ドメインの位置情報を手動でPfamのENV HXB2のページからとってきて，toml形式とかで保存しておく。

#memo
#irisでやる設計
#まずはファイル名等ある程度ハードコーディングで。その後アップデートしていく

#prerequisits
# python3 Anaconda 
# pip install biopython toml seaborn

#なぜかA,B,Cのみになってる
NAME_CORE=HIV-1-gM-A,B,C-noRs_env-aa
FASTA_FILE="../data/${NAME_CORE}.fasta"
HXB2_FILE="../data/HXB2_P04578.fasta"
DOMAIN_TOML_FILE="./toml/HIV-1_env_domains_HXB2.toml"
WD=`pwd`

if [[ ! -f ${FASTA_FILE} ]]; then
    echo "can't find ${FASTA_FILE}, exitting.."
    exit 1
fi
if [[ ! -f ${HXB2_FILE} ]]; then
    echo "can't find ${HXB2_FILE}, exitting.."
    exit 1
fi
if [[ ! -f ${DOMAIN_TOML_FILE} ]]; then
    echo "can't find ${DOMAIN_TOML_FILE}, exitting.."
    exit 1
fi

###########################
# extract domain sequences
###########################
#cat "${HXB2_FILE}" "${FASTA_FILE}" > "../data/with-HXB2_${NAME_CORE}.fasta"
#
#### MAFFT
#MAFFT_WORK_IRIS=/work/shohei1029/16spp/analysis/mafft
MAFFT_IN_FILE_NAME="with-HXB2_${NAME_CORE}.fasta"
MAFFT_OUT_FILE_NAME="mafft-linsi_${MAFFT_IN_FILE_NAME}"
#
#cp ../data/${MAFFT_IN_FILE_NAME} ${MAFFT_WORK_IRIS}
#cd $MAFFT_WORK_IRIS
#
#mafft-linsi --thread 128 $MAFFT_IN_FILE_NAME > $MAFFT_OUT_FILE_NAME
#
#cp $MAFFT_OUT_FILE_NAME ${WD}/../data/
#cd $WD
#
#jobdone mafft $MAFFT_IN_FILE_NAME @iris
####
#
#bioawk -cfastx '{if ($name !~ /HV1H2/){print ">"$name"\n"$seq}}' ../data/${MAFFT_OUT_FILE_NAME} > ../data/removed-HXB2_${MAFFT_OUT_FILE_NAME}
#↑先に実行してたはず

python extractDomain_getPosinToml_inAlignedFasta_createSecFasta.py -t ${DOMAIN_TOML_FILE} -f ../data/removed-HXB2_${MAFFT_OUT_FILE_NAME} > ../data/domains-HXB2-Pfam_${NAME_CORE}.fasta

mkdir -p ../data/domains_${NAME_CORE}

DOMAIN_LIST=(`python get_keys_inToml.py ${DOMAIN_TOML_FILE}`)
echo $DOMAIN_LIST

for DOMAIN in ${DOMAIN_LIST}
do
    bioawk -cfastx "{if (\$name ~ /-${DOMAIN}/){print \">\"\$name\"\n\"\$seq}}" ../data/domains-HXB2-Pfam_${NAME_CORE}.fasta > ../data/domains_${NAME_CORE}/${DOMAIN}_domains-HXB2-Pfam_${NAME_CORE}.fasta
    nohup mafft-linsi --thread 32 ../data/domains_${NAME_CORE}/${DOMAIN}_domains-HXB2-Pfam_${NAME_CORE}.fasta > ../data/domains_${NAME_CORE}/mafft_linsi_${DOMAIN}_domains-HXB2-Pfam_${NAME_CORE}.fasta &

    nohup zsh ./blastp_makesim_args.sh ${NAME_CORE} ${DOMAIN} &
done

jobdone "shell script ${NAME_CORE} (BLAST may be still running.)"
