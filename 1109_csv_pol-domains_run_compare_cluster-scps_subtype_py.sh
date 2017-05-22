#!/bin/zsh
#set -euo pipefail

#1109. 全simへ
#列を各領域として，CSVで表示

#echo -n で改行なし

#header
echo "seq identity threshold,RVP,RVT_1,RVT_thumb,RVT_connect,RNase_H,Integrase_Zn,rve,IN_DBD_C"
for i in {0..99}
do
    echo -n ${i}
    for DOM in RVP RVT_1 RVT_thumb RVT_connect RNase_H Integrase_Zn rve IN_DBD_C
    do
        echo -n ","
        tmp=`python compare_cluster-scps_subtype.py ~/ttckfs/16spp/analysis/clusterx/pol-domains-HXB2-Pfam/scps-c3_${DOM}_pi${i}_blastp_1e-5_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.txt`
        echo -n ${tmp} | tr '\n' ' '
    done
    echo ''
done


