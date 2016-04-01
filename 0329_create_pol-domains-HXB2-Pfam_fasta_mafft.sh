#!/bin/env zsh
set -euo pipefail

mkdir -p ../data/pol-domains-HXB2-Pfam
cd ../data/pol-domains-HXB2-Pfam

bioawk -cfastx '{if ($name ~ /-RVP/){print ">"$name"\n"$seq}}' ../domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > RVP_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 RVP_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta >  mafft-linsi_RVP_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &
bioawk -cfastx '{if ($name ~ /-RVT_1/){print ">"$name"\n"$seq}}' ../domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > RVT_1_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 RVT_1_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta >  mafft-linsi_RVT_1_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &
bioawk -cfastx '{if ($name ~ /-RVT_thumb/){print ">"$name"\n"$seq}}' ../domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > RVT_thumb_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 RVT_thumb_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta >  mafft-linsi_RVT_thumb_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &
bioawk -cfastx '{if ($name ~ /-RVT_connect/){print ">"$name"\n"$seq}}' ../domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > RVT_connect_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 RVT_connect_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta >  mafft-linsi_RVT_connect_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &
bioawk -cfastx '{if ($name ~ /-RNase_H/){print ">"$name"\n"$seq}}' ../domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > RNase_H_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 RNase_H_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta >  mafft-linsi_RNase_H_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &
bioawk -cfastx '{if ($name ~ /-Integrase_Zn/){print ">"$name"\n"$seq}}' ../domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > Integrase_Zn_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 Integrase_Zn_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta >  mafft-linsi_Integrase_Zn_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &
bioawk -cfastx '{if ($name ~ /-rve/){print ">"$name"\n"$seq}}' ../domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > rve_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 rve_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta >  mafft-linsi_rve_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &
bioawk -cfastx '{if ($name ~ /-IN_DBD_C/){print ">"$name"\n"$seq}}' ../domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta > IN_DBD_C_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta && mafft-linsi --thread 32 IN_DBD_C_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta >  mafft-linsi_IN_DBD_C_domains-HXB2-Pfam_A,B,C_HIV-1-gM-noRs_pol-aa_v3.fasta &





