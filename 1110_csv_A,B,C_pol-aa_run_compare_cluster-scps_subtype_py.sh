#!/bin/zsh
#set -euo pipefail

#localで実行。ttckfsをマウントした状態で。


#1109. 全simへ
#列を各領域として，CSVで表示
#1110 pol-aaへ

#echo -n で改行なし

#header
echo 'seq identity threshold,"HIV-1-gM-A,B,C-noRs_pol-aa"'
DOM="A,B,C"
for i in {50..75..5} {72..78..2} {80..99}
do
    echo -n ${i}
#    for DOM in RVP RVT_1 RVT_thumb RVT_connect RNase_H Integrase_Zn rve IN_DBD_C
#    do
        echo -n ","
        tmp=`python compare_cluster-scps_subtype.py ~/ttckfs/16spp/analysis/clusterx/HIV-1-gM-A,B,C-noRs_pol-aa/scps-c3_${DOM}_pi${i}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt`
        echo -n ${tmp} | tr '\n' ' '
#    done
    echo ''
done
