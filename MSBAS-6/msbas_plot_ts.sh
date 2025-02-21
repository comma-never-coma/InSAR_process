#! /bin/bash

[[ $# != 1  || ! -f ${1} ]] && {
    echo "provide file to plot, exitinig..."
    exit 1
}

nc=$(awk 'NR == 1 { print NF; exit; }' ${1})
outname=$(echo ${1} | awk '{gsub(/.txt/,"")}; 1')
read -r xmin xmax <<< $(awk 'BEGIN {min = 10000; max = 0} {if ($2<min) min=$2; if ($2>max) max=$2} END {print min-0.025*(max-min), max+0.025*(max-min)}' ${1})
error_bars=1

#set title '${1}' offset 0,0
#set bars small 

cat << EOF > ${1}.gnuplot
set xlabel "Time, year"
set xlabel offset 0,0.5
set ylabel offset 0.5,0
set ylabel "Displacement, m"
set autoscale
set xrange [${xmin}:${xmax}]
set size 1, 1
set key reverse Left bottom left samplen 2
set output "${outname}.pdf" 
set terminal pdfcairo enhanced crop font "Helvetica,16"
EOF

if [[ ${nc} -eq 4 ]]; then
if [[ ${error_bars} == 1 ]]; then
cat <<EOF >> ${1}.gnuplot
plot \
"${1}" using 2:3:4 with yerr title "Line-of-Sight" ps 0.2 pt 7 lw 1  lc rgb '#5DA5DA',\
"${1}" using 2:3 with lines notitle lw 1  lt 1 lc rgb '#5DA5DA'
EOF
else
cat << EOF >> ${1}.gnuplot
plot \
"${1}" using 2:3 with points title "Line-of-Sight" ps 0.4 pt 8 lc rgb '#60BD68'
EOF
fi    
fi 

if [[ ${nc} -eq 6 ]]; then
if [[ ${error_bars} == 1 ]]; then
cat <<EOF >> ${1}.gnuplot
plot \
"${1}" using 2:3:4 with yerr title "East" ps 0.2 pt 7 lw 1  lc rgb '#5DA5DA',\
"${1}" using 2:5:6 with yerr title "Vertical" ps 0.2 pt 7 lw 1  lc rgb '#B276B2',\
"${1}" using 2:3 with lines notitle lw 1  lt 1 lc rgb '#5DA5DA',\
"${1}" using 2:5 with lines notitle lw 1  lt 1 lc rgb '#B276B2'
EOF
else
cat << EOF >> ${1}.gnuplot
plot \
"${1}" using 2:3 with points title "East" ps 0.4 pt 4 lc rgb '#5DA5DA',\
"${1}" using 2:5 with points title "Vertical" ps 0.4 pt 12 lc rgb '#B276B2'
EOF
fi    
fi 

if [[ ${nc} -eq 8 ]]; then
if [[ ${error_bars} == 1 ]]; then
cat <<EOF >> ${1}.gnuplot
plot \
"${1}" using 2:3:4 with yerr title "North" ps 0.2 pt 7 lw 1  lc rgb '#60BD68',\
"${1}" using 2:5:6 with yerr title "East" ps 0.2 pt 7 lw 1  lc rgb '#5DA5DA',\
"${1}" using 2:7:8 with yerr title "Vertical" ps 0.2 pt 7 lw 1  lc rgb '#B276B2',\
"${1}" using 2:3 with lines notitle lw 1  lt 1 lc rgb '#60BD68',\
"${1}" using 2:5 with lines notitle lw 1  lt 1 lc rgb '#5DA5DA',\
"${1}" using 2:7 with lines notitle lw 1  lt 1 lc rgb '#B276B2'
EOF
else
cat << EOF >> ${1}.gnuplot
plot \
"${1}" using 2:3 with points title "North" ps 0.4 pt 8 lc rgb '#60BD68',\
"${1}" using 2:5 with points title "East" ps 0.4 pt 4 lc rgb '#5DA5DA',\
"${1}" using 2:7 with points title "Vertical" ps 0.4 pt 12 lc rgb '#B276B2'
EOF
fi    
fi 

if [[ ${nc} -eq 10 ]]; then
if [[ ${error_bars} == 1 ]]; then
cat <<EOF >> ${1}.gnuplot
plot \
"${1}" using 2:3:4 with yerr title "North" ps 0.2 pt 7 lw 1 lc rgb '#60BD68',\
"${1}" using 2:5:6 with yerr title "East" ps 0.2 pt 7 lw 1 lc rgb '#5DA5DA',\
"${1}" using 2:7:8 with yerr title "Vertical SPF" ps 0.2 pt 7 lw 1 lc rgb '#B276B2',\
"${1}" using 2:9:10 with yerr title "Vertical nSPF" ps 0.2 pt 7 lw 1 lc rgb '#F15854',\
"${1}" using 2:3 with lines notitle lw 1  lt 1 lc rgb '#60BD68',\
"${1}" using 2:5 with lines notitle lw 1  lt 1 lc rgb '#5DA5DA',\
"${1}" using 2:7 with lines notitle lw 1  lt 1 lc rgb '#B276B2',\
"${1}" using 2:9 with lines notitle lw 1  lt 1 lc rgb '#F15854'
EOF
else
cat << EOF >> ${1}.gnuplot
plot \
"${1}" using 2:3 with points title "North" ps 0.4 pt 8 lc rgb '#60BD68',\
"${1}" using 2:5 with points title "East" ps 0.4 pt 4 lc rgb '#5DA5DA',\
"${1}" using 2:7 with points title "Vertical SPF" ps 0.4 pt 12 lc rgb '#B276B2',\
"${1}" using 2:9 with points title "Vertical nSPF" ps 0.4 pt 6 lc rgb '#F15854'
EOF
fi    
fi 

gnuplot ${1}.gnuplot
rm ${1}.gnuplot

#    4D4D4D (gray)
#    5DA5DA (blue)
#    FAA43A (orange)
#    60BD68 (green)
#    F17CB0 (pink)
#    B2912F (brown)
#    B276B2 (purple)
#    DECF3F (yellow)
#    F15854 (red)
