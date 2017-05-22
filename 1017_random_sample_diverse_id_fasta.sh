#!/bin/env zsh
set -euo pipefail

# Created on 2016.10.17 by Shohei Nagata.

# ランダムに配列をサンプリングするスクリプトを実行して，
# 実行結果のfastaのID名の4文字分の重複を調べて，指定数以上のID名多様性が得られるまでサンプリングを繰り返す.


#TARGET_CORE="mafft-linsi_RVT_connect_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3"
#echo ${TARGET_CORE}
#while true
#do
#    IDVAR=`cat ../data/pol-domains-HXB2-Pfam/${TARGET_CORE}.fasta | python random_sample_fasta_subtype.py | tee >(python fastaID_diversity.py) > ../data/pol-domains-HXB2-Pfam/10seqseach/trials/10seqs_${TARGET_CORE}.t.fasta`
#    if [ ${IDVAR} = "30" ]; then
#        break
#    fi
#done



TARGET_CORE="mafft-linsi_RNase_H_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3"
echo ${TARGET_CORE}
while true
do
    IDVAR=`cat ../data/pol-domains-HXB2-Pfam/${TARGET_CORE}.fasta | python random_sample_fasta_subtype.py | tee >(python fastaID_diversity.py) > ../data/pol-domains-HXB2-Pfam/10seqseach/trials/10seqs_${TARGET_CORE}.t.fasta`
    echo ${IDVAR}
    if [ ${IDVAR} = "30" ]; then
        break
    fi
done

