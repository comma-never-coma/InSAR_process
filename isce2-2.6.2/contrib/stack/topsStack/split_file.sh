#!/bin/bash
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo " `basename $0`: split run_files "
echo   "17-07-2023 by hao.dou"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
if [ $# -ne 1 ]; then
        echo " split_file.sh [parallel_num] "
        echo " parallel_num (numbers of cores) "
        exit 1 
fi

parallel_num=$1
files=("run_02_unpack_secondary_slc" "run_05_overlap_geo2rdr"
"run_06_overlap_resample" "run_07_pairs_misreg" "run_09_fullBurst_geo2rdr"
"run_10_fullBurst_resample" "run_13_generate_burst_igram"
"run_14_merge_burst_igram" "run_15_filter_coherence" "run_16_unwrap")

for input_file in "${files[@]}"; do
    total_lines=$(wc -l < ${input_file})
    lines_per_file=$(( ($total_lines + $parallel_num - 1) / $parallel_num ))
    awk -v input_file="${input_file}" -v lines_per_file=${lines_per_file} 'NR%lines_per_file == 1 {file_num++; file=sprintf("%s_%01d", input_file, file_num);} {print > file}' ${input_file}
done
