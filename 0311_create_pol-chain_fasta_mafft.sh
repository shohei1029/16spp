#!/bin/env zsh
set -euo pipefail

cd ../data
#bioawk -cfastx '{if ($name ~ /-Integrase/){print ">"$name"\n"$seq}}' pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > Integrase_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 Integrase_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta >  mafft-linsi_Integrase_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &

bioawk -cfastx '{if ($name ~ /-p51_RT/){print ">"$name"\n"$seq}}' pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > p51_RT_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 p51_RT_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta >  mafft-linsi_p51_RT_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &

bioawk -cfastx '{if ($name ~ /-p15/){print ">"$name"\n"$seq}}' pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > p15_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 p15_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > mafft-linsi_p15_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &

bioawk -cfastx '{if ($name ~ /-Protease/){print ">"$name"\n"$seq}}' pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > Protease_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 Protease_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > mafft-linsi_Protease_pol-chains_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &


