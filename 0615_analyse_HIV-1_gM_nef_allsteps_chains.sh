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
#とりあえず，サブタイプをA, B, Cに絞らずやってみる
#stderrをログへ出力する形で実行したいね
#そのうち，swissprotのtxtファイルを読みこませるだけで，どの領域とってきて〜ってのから自動でできるようにしたいかも。でもdomain使うかchainデータ使うかとか目で見ていきたさがある。

#☆手動で設定すること
#・TARGET_NAMEとPROT_LENGTH
#・domain/chain等の情報をtomlに書く
#・HXB2での目標タンパク質の配列をとってきて../data/内に置いておく

TARGET_NAME=nef
PROT_LENGTH=206
NAME_CORE=HIV-1-gM_${TARGET_NAME}-aa
FASTA_FILE_NAME=${NAME_CORE}.fasta
FASTA_FILE="../data/${FASTA_FILE_NAME}" 
WD=`pwd`

TARGET_REGION_TYPE=chains #e.g. domains, chains
HXB2_FILE="../data/HXB2_P04601.fasta" #nef
DOMAIN_TOML_FILE="./toml/HIV-1_${TARGET_NAME}_${TARGET_REGION_TYPE}_HXB2.toml"

if [[ ! -f ${HXB2_FILE} ]]; then
    echo "can't find ${HXB2_FILE}, exitting.."
    exit 1
fi
if [[ ! -f ${DOMAIN_TOML_FILE} ]]; then
    echo "can't find ${DOMAIN_TOML_FILE}, exitting.."
    exit 1
fi

########################
# make data sets #
########################
#-t:target gene/protein name , -m s: specify subtype (remove recombinants) mode , -s: specify subtype (sep=,)| -n: remove gene name in fasta header, criteria default: 90%, refseqlength: isolate HXB2
python getseq_ingbk_withbginfo_tofasta.py -t ${TARGET_NAME} -m o -o HIV-1 -e N,O,P | python filter_fasta_by_seq-length.py -n True -r ${PROT_LENGTH} > $FASTA_FILE
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
for pi in 0 {5..6}{0,5} {7..9}{0..9}
do
    cat ../analysis/blast/mbs_${BLAST_OUT_FILE_NAME} | python ./v5_blast7_tosims.py -p $pi > ../analysis/sim/${NAME_CORE}/pi${pi}_${BLAST_OUT_FILE_NAME}
done

jobdone create sif/sim file ${NAME_CORE} @iris, now you can visualise your network :D

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
nohup scp -r ../analysis/sim/${NAME_CORE} shohei@133.27.17.109:/media/shohei1029/16spp/analysis/sim &




# --領域ごとの配列類似性ネットワーク構築スクリプト，連結--
# もち修正含む。(ファイル転送をnohupにしたり)
#prerequisits
# python3 Anaconda 
# pip install biopython toml seaborn

#なぜかA,B,Cのみになってる
#NAME_CORE=HIV-1-gM-A,B,C-noRs_env-aa
#FASTA_FILE="../data/${NAME_CORE}.fasta"

if [[ ! -f ${FASTA_FILE} ]]; then
    echo "can't find ${FASTA_FILE}, exitting.."
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

cp ../data/${MAFFT_IN_FILE_NAME} ${MAFFT_WORK_IRIS}
cd $MAFFT_WORK_IRIS

mafft-linsi --thread 128 $MAFFT_IN_FILE_NAME > $MAFFT_OUT_FILE_NAME

cp $MAFFT_OUT_FILE_NAME ${WD}/../data/
cd $WD

jobdone mafft $MAFFT_IN_FILE_NAME @iris
###

bioawk -cfastx '{if ($name !~ /HV1H2/){print ">"$name"\n"$seq}}' ../data/${MAFFT_OUT_FILE_NAME} > ../data/removed-HXB2_${MAFFT_OUT_FILE_NAME}
#↑先に実行してたはず

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

    nohup zsh ./blastp_makesim_args.sh ${NAME_CORE} ${DOMAIN} &
done

jobdone "shell script ${NAME_CORE} (BLAST may be still running.)"
