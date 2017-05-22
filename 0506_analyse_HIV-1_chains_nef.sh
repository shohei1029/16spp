#!/bin/env zsh
set -euo pipefail

# 2016.4.24 by Shohei Nagata
# ドメインごとの領域を対象にsimファイルの生成まで。
# ドメインはisolate HXB2を基準に。ドメインの位置情報を手動でPfamのENV HXB2のページからとってきて，toml形式とかで保存しておく。

# 2016.5.4
#domainsからchainsへ変更．
#HXB2を付けてのアライメント等は既に実行されているはずだからコメントアウト

# 2016.5.9
#vifとかでやったときのshを参考に色々アップデート

#memo
#irisでやる設計
#まずはファイル名等ある程度ハードコーディングで。その後アップデートしていく

#prerequisits
# python3 Anaconda 
# pip install biopython toml seaborn

#☆手動で設定すること
#・TARGET_NAMEとPROT_LENGTH
#・domain/chain等の情報をtomlに書く
#・HXB2での目標タンパク質の配列をとってきて../data/内に置いておく

TARGET_NAME=nef
NAME_CORE=HIV-1-gM-A,B,C-noRs_${TARGET_NAME}-aa
FASTA_FILE_NAME=${NAME_CORE}.fasta
FASTA_FILE="../data/${FASTA_FILE_NAME}" 
WD=`pwd`

TARGET_REGION_TYPE=chains #e.g. domains, chains
HXB2_FILE="../data/HXB2_P04601.fasta"
DOMAIN_TOML_FILE="./toml/HIV-1_${TARGET_NAME}_${TARGET_REGION_TYPE}_HXB2.toml"

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
cat "${HXB2_FILE}" "${FASTA_FILE}" > "../data/with-HXB2_${NAME_CORE}.fasta"

### MAFFT
MAFFT_WORK_IRIS=/work/shohei1029/16spp/analysis/mafft
MAFFT_IN_FILE_NAME="with-HXB2_${NAME_CORE}.fasta"
MAFFT_OUT_FILE_NAME="mafft-linsi_${MAFFT_IN_FILE_NAME}"

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
##↑先に実行してたはず

python extractDomain_getPosinToml_inAlignedFasta_createSecFasta.py -t ${DOMAIN_TOML_FILE} -f ../data/removed-HXB2_${MAFFT_OUT_FILE_NAME} > ../data/${TARGET_REGION_TYPE}-HXB2-Pfam_${NAME_CORE}.fasta

mkdir -p ../data/${TARGET_REGION_TYPE}_${NAME_CORE}

DOMAIN_LIST=(`python get_keys_inToml.py ${DOMAIN_TOML_FILE}`)
echo $DOMAIN_LIST

for DOMAIN in ${DOMAIN_LIST}
do
    bioawk -cfastx "{if (\$name ~ /-${DOMAIN}/){print \">\"\$name\"\n\"\$seq}}" ../data/${TARGET_REGION_TYPE}-HXB2-Pfam_${NAME_CORE}.fasta > ../data/${TARGET_REGION_TYPE}_${NAME_CORE}/${DOMAIN}_${TARGET_REGION_TYPE}-HXB2-Pfam_${NAME_CORE}.fasta
    nohup mafft-linsi --thread 32 ../data/${TARGET_REGION_TYPE}_${NAME_CORE}/${DOMAIN}_${TARGET_REGION_TYPE}-HXB2-Pfam_${NAME_CORE}.fasta > ../data/${TARGET_REGION_TYPE}_${NAME_CORE}/mafft_linsi_${DOMAIN}_${TARGET_REGION_TYPE}-HXB2-Pfam_${NAME_CORE}.fasta &

    nohup zsh ./blastp_makesim_args_v2.sh ${NAME_CORE} ${DOMAIN} ${TARGET_REGION_TYPE} &
done

jobdone "shell script ${NAME_CORE} (BLAST may be still running.)"
