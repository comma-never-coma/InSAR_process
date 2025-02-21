#!/bin/bash

#print help information if no parameters are given
function usage() {
	cat <<END
*** $(basename $0)
*** Copyright 2018, NRCan

usage: $(basename $0) -[c:fil:s] <file_tif>

input parameters (options must preceed arguments):
	-c [value]              use provided color palette  
	-f                      retain full resolution
	-i                      clipt min/max to [mean-5std,mean+5std]
	-l [value1:value2] 		clip min/max to [value1,value2], value1<value2
	-s                      use symmetric color palette  
	file_tif                (input) *tif file (in radar coordinates or geocoded, individual color palettes are used for files containing keywords: wrp, dem, cc, mli)
END
}

fullResolution=0
colorPalette=0
symmetricPalette=0
intellegentPalette=0
definedPalette=0
while getopts "c:fil:s" OPTION; do
	case ${OPTION} in
	c)
		colorPalette=${OPTARG,,}
		;;
	f)
		fullResolution=1
		;;
	i)
		intelligentPalette=1
		;;
	l)
		IFS=: read -r -a definedPalette <<< "${OPTARG}"
		;;
	s)
		symmetricPalette=1
		;;
	?)
		duapMsg "warning in $(basename ${0}): unexpected option was ignored"
		;;
	esac
done
shift $((OPTIND - 1))

#validate input paramets
[[ $# -eq 0 || ! -s ${1} ]] && {
	duapMsg "error in $(basename $0): input parameter(s) are missing or incorrect, exiting..."
	usage
	exit 1
}

return_dir=$(pwd -L)
work_dir=$(dirname ${1})

cd ${work_dir}
infile=$(basename ${1}).grd
outfile=$(echo ${infile} | awk '{gsub(/.tif.grd/,"")}; 1').ps
cpt=$(echo ${infile} | awk '{gsub(/.tif.grd/,"")}; 1').cpt
[[ -s ${cpt} ]] && rm ${cpt}

# convert zeros to NaNs, tif to grd
gmt grdmath $(basename ${1}) 0 NAN = ${infile}.tmp

# need to know original (before resampling) xinc, yinc, xlength, and ylength
read -r xmin xmax xinc xlength <<<$(gmt grdinfo ${infile}.tmp | grep x_min: | awk '{print $3,$5,$7,$11}')
read -r ymin ymax yinc ylength <<<$(gmt grdinfo ${infile}.tmp | grep y_min: | awk '{print $3,$5,$7,$11}')

# if file is bigger than 2000x2000, resample to corser-resolution using nearest-neigbour interpolation unless full resolution is required
if [[ $(echo ${xlength} ${ylength} | awk '{flag=0; if ($1>2000 && $2>2000) {flag=1}; print flag;}') -eq 1 && ${fullResolution} -ne 1 ]]; then

	gmt grdsample -nn -I2000+/2000+ ${infile}.tmp -G${infile}
	rm ${infile}.tmp
else
	mv ${infile}.tmp ${infile}
fi

# define color palette based on file extension
if [[ ${colorPalette} != "0" ]]; then
	echo "applying provided palette"
elif [[ ${infile} == *"wrp"* ]]; then
	colorPalette="cyclic"
elif [[ ${infile} == *"dem"* ]]; then
	colorPalette="topo"
elif [[ ${infile} == *"cc"* ]]; then
	colorPalette="gray"
elif [[ ${infile} == *"mli"* ]]; then
	gmt grdmath -V ${infile} LOG10 10 MUL = ${infile}.tmp
	mv ${infile}.tmp ${infile}
	colorPalette="gray"
else
	colorPalette="polar"
fi

# define palette limits
if [[ ${intellegentPalette} != "0" ]]; then
	read -r zmean zstd <<<$(gmt grdinfo -L2 ${infile} | grep stdev: | awk '{print $3,$5}')
 	read -r zmin zmax <<<$(echo ${zmean} ${zstd} 5 | awk '{print $1-$2*$3,$1+$2*$3}')
elif [[ ${definedPalette} != "0" ]]; then
	read -r zmin zmax <<<$(echo ${definedPalette[0]} ${definedPalette[1]})
else
	read -r zmin zmax <<<$(gmt grdinfo ${infile} | grep z_min: | awk '{print $3,$5}')
fi

if [[ ${symmetricPalette} != "0" ]]; then
	read -r zmin zmax <<<$(echo ${zmin} ${zmax} | awk 'function abs(v) {return v < 0 ? -v : v} {if (abs($1)>abs($2)) print -abs($1),abs($1); else print -abs($2),abs($2);}')
fi

gmt makecpt -T${zmin}/${zmax}/256+ -C${colorPalette:-polar} -Z -D >${cpt}

range=${xmin}/${xmax}/${ymin}/${ymax}
pwidth=19 #cm

if [[ $(echo ${xinc} ${yinc} | awk '{flag=0; if ($1==1 && $2==1) {flag=1}; print flag;}') -eq 1 ]]; then

	# XY grid
	pscale=$(echo ${xmax} ${ymax} ${pwidth} | awk '{print $3*$2/$1}')
	projection=X${pwidth}c/${pscale}c
	gmt psbasemap -B+n -J$projection -R$range -X0.2c -Y0.2c -P -V -K --PROJ_LENGTH_UNIT=c --PS_MEDIA=A2 --MAP_FRAME_TYPE=plain --FORMAT_GEO_MAP=ddd:mm --MAP_FRAME_AXES=WeSn >${outfile} 
	gmt grdimage ${infile} -J$projection -R$range -O -V -P -C${cpt} -K >>${outfile}
	gmt psscale -D3.5c/1.5c/6c/0.5ch -C${cpt} -P -O -Ba >>${outfile}
	gmt psconvert ${outfile} -Tf -A

else

	# geographic grid, configure scale
	xmean=$(echo ${xmin} ${xmax} | awk '{print ($1+$2)/2}')
	ymean=$(echo ${ymin} ${ymax} | awk '{print ($1+$2)/2}')

	# Mercator projection
	#projection=M${xmin}/${ymin}/${pwidth}c

	# Lambert Conformal Conic projection - better for plotting large Subarctic areas
	projection=L${xmean}/${ymean}/${ymin}/${ymax}/${pwidth}c

	scale=2560 #scale initial value 2560 km, bisect until it is less than 1/3 of the map, this will avoid non-integer values
	map_width_third=$(echo ${xmin} ${ymin} | gmt mapproject -J$projection -R$range -G${xmax}/${ymin}/k | awk '{print $3/3.0}')
	while [[ $(echo ${scale} ${map_width_third} | awk '{if ($1 > $2) print 1; else print 0;}') -eq "1" ]]; do scale=$(echo ${scale} | awk '{print int($1/2.0)}'); done

	gmt psbasemap -Bag -J$projection -R$range -X2.5c -Y2c -P -V -K --PROJ_LENGTH_UNIT=c --PS_MEDIA=A2 --MAP_FRAME_TYPE=plain --FORMAT_GEO_MAP=ddd:mm --MAP_FRAME_AXES=WeSn >${outfile}
	gmt pscoast -J$projection -R$range -Na -Df -O -V -W0.8,darkgray -S225/225/255 -G255/255/255 -P -K >>${outfile}
	gmt grdimage ${infile} -J$projection -R$range -O -V -P -C${cpt} -K -n+c -Q >>${outfile}
	gmt psbasemap -J$projection -R$range -O -P -V -K -Lx15c/1c+f+c${ymean}+w${scale}+l >>${outfile}

	# plot location of reference
	[[ -s ref_ll.region ]] && {
		gmt psxy -J$projection -R$range  -O -K -L -W1.5,black ref_ll.region >> ${outfile}
	}

	[[ -s ref_ll.point ]] && {
    	gmt pstext ref_ll.point -R$region -J$projection -O -K -P -D0.1c/0.1c -V -F+f16p,Helvetica,black+j+a >>${outfile}
	}

	[[ -s points.txt ]] && {
    	gmt psxy points.txt -R$region  -J$projection -W1,black -Sc0.15c -O -V -K >>${outfile} 
    	gmt pstext points.txt -R$region -J$projection -O -K -P -D0.1c/0.1c -V -F+f16p,Helvetica,black+j+a >>${outfile}
	}

	# prints file name
	#echo $1 | gmt pstext -J$projection  -R$range -F+cTL -Gwhite -W0.8,darkgray -O -P -K >> ${outfile}

	gmt psscale -D3.5c/1.5c/6c/0.5ch -C${cpt} -P -O -Ba >>${outfile} #-Ba+lm - add unit m
	gmt psconvert ${outfile} -Tf -A

	# plot kml
	gmt grdimage ${infile} -JX10d -R$range --MAP_FRAME_TYPE=inside -C${cpt} -Q >${outfile}
	#gmt psscale -D3.5c/1.5c/6c/0.5ch -C${cpt}  -P -O -Ba >> ${outfile}
	gmt psconvert ${outfile} -TG -W+k+t"$(basename ${infile} .grd)"+n"$(basename ${infile} .grd)"+l256/-1

	png=$(echo ${infile} | awk '{gsub(/.tif.grd/,"")}; 1').png
	kml=$(echo ${infile} | awk '{gsub(/.tif.grd/,"")}; 1').kml
	kmz=$(echo ${infile} | awk '{gsub(/.tif.grd/,"")}; 1').kmz

	zip -j ${kmz} ${kml} ${png}
	rm ${kml} ${png}
fi

rm ${cpt} ${infile} ${outfile} #${gmtconf}
[[ -s gmt.history ]] && rm gmt.history

cd ${return_dir}
