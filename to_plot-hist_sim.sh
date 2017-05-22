#!/bin/env zsh
set -euo pipefail

# 2016.6.6 by Shohei Nagata
# find . -name "*.txt" | xargs -n 1 sh to_plot_sim.sh

echo "$@"

input_path=$@
input_filename=${input_path##*/}
input_filename_without_extension=${input_filename%.*}
input_dir_path=${input_path%/*}

#echo ${input_filename}
#echo ${input_filename_without_extension}
#echo ${input_dir_path}

cat ${input_path} | cut -f 3 -d " " | python ./plot_hist.py -o ${input_dir_path}/${input_filename_without_extension}.png --title ${input_filename_without_extension} --xlim 0,1 --xlabel "sim" --ylabel "counts" 
#cat ${input_path} | cut -f 3 -d " " | python ./plot_hist.py -o ${input_dir_path}/${input_filename_without_extension}.png --title ${input_filename_without_extension} --xlim 0,1 --ylim 0,5000 --xlabel "sim" --ylabel "counts" 
