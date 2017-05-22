#!/bin/zsh
#set -euo pipefail

# pi 87 91 93 97 99

perl -lane 'print "$F[0] $F[1] $F[2]" if ($F[2] >= 0.87)' ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi0_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt > ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi87_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt
PI=87
clusterx -c 3 -o ../analysis/clusterx/HIV-1-gM-noRs_pol-aa/scps-c3_pi${PI}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi${PI}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt &

perl -lane 'print "$F[0] $F[1] $F[2]" if ($F[2] >= 0.91)' ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi0_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt > ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi91_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt
PI=91
clusterx -c 3 -o ../analysis/clusterx/HIV-1-gM-noRs_pol-aa/scps-c3_pi${PI}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi${PI}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt &

perl -lane 'print "$F[0] $F[1] $F[2]" if ($F[2] >= 0.93)' ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi0_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt > ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi93_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt
PI=93
clusterx -c 3 -o ../analysis/clusterx/HIV-1-gM-noRs_pol-aa/scps-c3_pi${PI}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi${PI}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt &

perl -lane 'print "$F[0] $F[1] $F[2]" if ($F[2] >= 0.97)' ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi0_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt > ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi97_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt
PI=97
clusterx -c 3 -o ../analysis/clusterx/HIV-1-gM-noRs_pol-aa/scps-c3_pi${PI}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi${PI}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt &

perl -lane 'print "$F[0] $F[1] $F[2]" if ($F[2] >= 0.99)' ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi0_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt > ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi99_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt
PI=99
clusterx -c 3 -o ../analysis/clusterx/HIV-1-gM-noRs_pol-aa/scps-c3_pi${PI}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt ../analysis/sim/HIV-1-gM-noRs_pol-aa/pi${PI}_blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.sim.txt &
