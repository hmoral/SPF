#!/bin/bash
Nin=$1
opt1=$2
optDif=$3
adapP=$4
seed=$5
collapseN=$6
outPref_in=$7
U=$8
adapU=$9

coord_file_in=coords_SLiM_GCF_000247815.1_FicAlb1.5_chr23.txt
rho_file_in=out_rho4Ne_win2e+05_col.hun.snpPair_runALL_chr23.txt
hs_in=1_hs_prod1_hmode1_Ne20K_step40K.txt
g=5000

echo "START at " $(date '+%d_%m_%Y_%H_%M_%S') "COMMAND:  " slim -t -d "adapU=${adapU}" -d "U=${U}" -d "g=${g}" -d "outPref_in='${outPref_in}'" -d "coord_file_in='${coord_file_in}'" -d "rho_file_in='${rho_file_in}'" -d "hs_in='${hs_in}'" -d "collapseN=${collapseN}" -d "Nin=${Nin}" -d "opt1=${opt1}" -d "optDif=${optDif}" -d "seed=${seed}" -d "adapP=${adapP}" SPF_loadStabSel_3_10_22.slim
slim -t -d "adapU=${adapU}" -d "U=${U}" -d "g=${g}" -d "outPref_in='${outPref_in}'" -d "coord_file_in='${coord_file_in}'" -d "rho_file_in='${rho_file_in}'" -d "hs_in='${hs_in}'" -d "collapseN=${collapseN}" -d "Nin=${Nin}" -d "opt1=${opt1}" -d "optDif=${optDif}" -d "seed=${seed}" -d "adapP=${adapP}" SPF_loadStabSel_3_10_22.slim
echo "END at " $(date '+%d_%m_%Y_%H_%M_%S') "COMMAND:  " slim -t -d "adapU=${adapU}" -d "U=${U}" -d "g=${g}" -d "outPref_in='${outPref_in}'" -d "coord_file_in='${coord_file_in}'" -d "rho_file_in='${rho_file_in}'" -d "hs_in='${hs_in}'" -d "collapseN=${collapseN}" -d "Nin=${Nin}" -d "opt1=${opt1}" -d "optDif=${optDif}" -d "seed=${seed}" -d "adapP=${adapP}" SPF_loadStabSel_3_10_22.slim
