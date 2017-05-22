#!/bin/env zsh
set -euo pipefail

# 2016.4.11 by Shohei Nagata
# env領域のアミノ酸配列抽出，総当りblast，sim作成，mafftによるアライメントまで
# v1.0.1
#  生成されるsimファイル名のバグを修正
# 4.26
# vif 用へ。ドメイン等なしで。
# その後
# tatとか。他のタンパク質でも使えるぜ。
# --連結--
# ドメインごとの領域を対象にsimファイルの生成まで。
# ドメインはisolate HXB2を基準に。ドメインの位置情報を手動でPfamのENV HXB2のページからとってきて，toml形式とかで保存しておく。


#memo
#irisでやる設計
#stderrをログへ出力する形で実行したいね
#そのうち，swissprotのtxtファイルを読みこませるだけで，どの領域とってきて〜ってのから自動でできるようにしたいかも。でもdomain使うかchainデータ使うかとか目で見ていきたさがある。

#☆手動で設定すること
#・TARGET_NAMEとPROT_LENGTH
#・domain/chain等の情報をtomlに書く
#・HXB2での目標タンパク質の配列をとってきて../data/内に置いておく.file名も設定する。
#・TARGET_REGION_TYPEの設定

TARGET_NAME=pol
PROT_LENGTH=1003 #gag-polの1435からgagの500を引いただけ。pol単体ではswiss-protに記載なし。 →いや，だめだ。前みたときの1003をもとにしよう。
NAME_CORE=A,B,C_HIV-1-gM-noRs_${TARGET_NAME}-aa_v3
FASTA_FILE_NAME=${NAME_CORE}.fasta
FASTA_FILE="../data/${FASTA_FILE_NAME}" 
WD=`pwd`

TARGET_REGION_TYPE=non-domains-Pfam #e.g. domains, chains
HXB2_FILE="../data/HXB2_P04585.fasta"  #gag-pol
#HXB2_FILE="../data/HXB2_P04591.fasta"  #gag
#HXB2_FILE="../data/HXB2_P04578.fasta" #env
DOMAIN_TOML_FILE="./toml/HIV-1_${TARGET_NAME}_${TARGET_REGION_TYPE}_HXB2.toml"

#if [[ ! -f ${HXB2_FILE} ]]; then
#    echo "can't find ${HXB2_FILE}, exitting.."
#    exit 1
#fi
if [[ ! -f ${DOMAIN_TOML_FILE} ]]; then
    echo "can't find ${DOMAIN_TOML_FILE}, exitting.."
    exit 1
fi

#なぜかA,B,Cのみになってる
#NAME_CORE=HIV-1-gM-A,B,C-noRs_env-aa
#FASTA_FILE="../data/${NAME_CORE}.fasta"

#if [[ ! -f ${FASTA_FILE} ]]; then
#    echo "can't find ${FASTA_FILE}, exitting.."
#    exit 1
#fi

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
##↑先に実行してたはず

#このスクリプト(0517_gag)より追加。swiss-protを基にしたtomlに書かれた位置を，アライメント後fastaの1番目の配列(基本的にHXB2にしてる)上の位置へ変換して(ギャップ分変換)，新しくtomlを吐き出すスクリプト.
DOMAIN_TOML_FILE_ALIGNED="./toml/aligned_HIV-1_${TARGET_NAME}_${TARGET_REGION_TYPE}_HXB2.toml"
python convert_toml-seqpos_toPosWithGap_1stSeqInFastaAsStandard.py -t ${DOMAIN_TOML_FILE} -f ../data/${MAFFT_OUT_FILE_NAME} > ${DOMAIN_TOML_FILE_ALIGNED}
#

DOMAIN_TOML_FILE=${DOMAIN_TOML_FILE_ALIGNED}
python extractDomain_getPosinToml_inAlignedFasta_createSecFasta.py -t ${DOMAIN_TOML_FILE} -f ../data/removed-HXB2_${MAFFT_OUT_FILE_NAME} > ../data/${TARGET_REGION_TYPE}-HXB2-Pfam_${NAME_CORE}.fasta

mkdir -p ../data/${TARGET_REGION_TYPE}_${NAME_CORE}

DOMAIN_LIST=(`python get_keys_inToml.py ${DOMAIN_TOML_FILE}`)
echo $DOMAIN_LIST

for DOMAIN in ${DOMAIN_LIST}
do
    bioawk -cfastx "{if (\$name ~ /-${DOMAIN}/){print \">\"\$name\"\n\"\$seq}}" ../data/${TARGET_REGION_TYPE}-HXB2-Pfam_${NAME_CORE}.fasta > ../data/${TARGET_REGION_TYPE}_${NAME_CORE}/${DOMAIN}_${TARGET_REGION_TYPE}-HXB2-Pfam_${NAME_CORE}.fasta
    nohup mafft-linsi --thread 32 ../data/${TARGET_REGION_TYPE}_${NAME_CORE}/${DOMAIN}_${TARGET_REGION_TYPE}-HXB2-Pfam_${NAME_CORE}.fasta > ../data/${TARGET_REGION_TYPE}_${NAME_CORE}/mafft_linsi_${DOMAIN}_${TARGET_REGION_TYPE}-HXB2-Pfam_${NAME_CORE}.fasta &

    nohup zsh ./blastp_makesim_args_v2.sh ${NAME_CORE} ${DOMAIN} ${TARGET_REGION_TYPE} & #v2はドメイン以外にも対応
done

jobdone "shell script ${NAME_CORE} (BLAST may be still running.)"
