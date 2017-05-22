#!/bin/zsh
#set -euo pipefail

for DOM in RVP RVT_1 RVT_thumb RVT_connect RNase_H Integrase_Zn rve IN_DBD_C
do
echo $DOM
    for i in {65..95}
    do
        python compare_cluster-scps_subtype.py ~/ttckfs/16spp/analysis/clusterx/pol-domains-HXB2-Pfam/scps-c3_${DOM}_pi${i}_blastp_1e-5_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.txt
    done
done


