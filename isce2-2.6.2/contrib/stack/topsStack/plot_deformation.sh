#! /bin/bash
##########################################################
# Script that:
# -plotting deformation maps
#       
#########################################################################################################

##########################################################################################################
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo " `basename $0`: plotting deformation  maps "
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
if [ $# -ne 6 ]; then
        echo " plot_timeseries <smallbaseline_directory> <output_directory> <cpt> [cpt_min] [cpt_max] [cpt_min_ts] [cpt_max_ts] "
        echo "      smallbaseline_directory: 和baselines、merged、reference、secondarys位于同级目录的时序处理工作目录 "
        echo "      output_directory: output .tiff files path "
        echo "      cpt_min: display minimum z value(mm/year), default: -1 (auto calculation)"
        echo "      cpt_max: display maximum z value(mm/year), default: 1 (auto calculation)"
        echo "      cpt_min_ts: display minimum z value(mm), default: -1 (auto calculation)"
        echo "      cpt_max_ts: display maximum z value(mm), default: 1 (auto calculation)"
        exit 1
fi

dir=$1
output_dir=$2

mkdir -vp $output_dir
cp $dir/pic/network.pdf $output_dir

cd $dir/geo && rm -rf PLOT
mkdir -vp PLOT && cd PLOT
save_gmt.py ../geo_avgSpatialCoh.h5 -o ./geo_avgSpatialCoh
save_gmt.py ../geo_temporalCoherence.h5 -o ./geo_temporalCoherence

gmt grdmath geo_avgSpatialCoh = geo_avgSpatialCoh.grd
gmt grdmath geo_temporalCoherence = geo_temporalCoherence.grd
gmt grd2xyz geo_avgSpatialCoh.grd > geo_avgSpatialCoh.xyz
sed -e '/NaN/d' geo_avgSpatialCoh.xyz > geo_avgSpatialCoh.xyz1
gmt grd2xyz geo_temporalCoherence.grd > geo_temporalCoherence.xyz
sed -e '/NaN/d' geo_temporalCoherence.xyz > geo_temporalCoherence.xyz1

lon1=`gdalinfo geo_avgSpatialCoh.grd | sed -n '/^Lower Left/p' | awk '{print $4}' | sed 's/,//g'`
lat1=`gdalinfo geo_avgSpatialCoh.grd | sed -n '/^Lower Left/p' | awk '{print $5}' | sed 's/)//g'`
lon2=`gdalinfo geo_avgSpatialCoh.grd | sed -n '/^Upper Right/p' | awk '{print $4}' | sed 's/,//g'`
lat2=`gdalinfo geo_avgSpatialCoh.grd | sed -n '/^Upper Right/p' | awk '{print $5}' | sed 's/)//g'`
lim=-R$lon1/$lon2/$lat1/$lat2

gmt gmtset FORMAT_GEO_MAP ddd:mm:ssF
gmt gmtset MAP_ANNOT_OFFSET_PRIMARY 0.07i MAP_ANNOT_OBLIQUE 34 MAP_FRAME_WIDTH 0.08i MAP_SCALE_HEIGHT 0.1i MAP_LABEL_OFFSET 0.08i MAP_TICK_LENGTH_PRIMARY 0.1i PS_MEDIA B3 MAP_FRAME_TYPE fancy COLOR_NAN white FONT 10p

gmt makecpt -Cjet -T0/1/0.01 -V -Z > g.cpt
gmt psbasemap $lim -JM15 -Bxaf -Byaf -BWeSn -P -K -V  > geo_avgSpatialCoh.ps
gmt psxy geo_avgSpatialCoh.xyz1 $lim -JM15 -Cg.cpt -V -Sc0.01i -K -O >> geo_avgSpatialCoh.ps
gmt psscale -Cg.cpt -Dx4.1i/0.35i+w1.65i/0.125i+h -I -O -Baf >> geo_avgSpatialCoh.ps
gmt psconvert geo_avgSpatialCoh.ps -Tg -A -E800 -W+g

gmt psbasemap $lim -JM15 -Bxaf -Byaf -BWeSn -P -K -V  > geo_temporalCoherence.ps
gmt psxy geo_temporalCoherence.xyz1 $lim -JM15 -Cg.cpt -V -Sc0.01i -K -O >> geo_temporalCoherence.ps
gmt psscale -Cg.cpt -Dx4.1i/0.35i+w1.65i/0.125i+h -I -O -Baf >> geo_temporalCoherence.ps
gmt psconvert geo_temporalCoherence.ps -Tg -A -E800 -W+g
mv *.tiff $output_dir && rm *


cd $dir/geo && rm -rf PLOT1
mkdir -vp PLOT1 && cd PLOT1
mask.py ../geo_velocity.h5 -m ../geo_maskTempCoh.h5 -o ./geo_velocity_msk.h5
save_gmt.py geo_velocity_msk.h5 -o geo_velocity

gmt grdmath geo_velocity 0 NAN 1000 MUL = geo_velocity.grd
gmt grd2xyz geo_velocity.grd > geo_velocity.xyz
sed -e '/NaN/d' geo_velocity.xyz > geo_velocity.xyz1

lon1=`gdalinfo geo_velocity.grd | sed -n '/^Lower Left/p' | awk '{print $4}' | sed 's/,//g'`
lat1=`gdalinfo geo_velocity.grd | sed -n '/^Lower Left/p' | awk '{print $5}' | sed 's/)//g'`
lon2=`gdalinfo geo_velocity.grd | sed -n '/^Upper Right/p' | awk '{print $4}' | sed 's/,//g'`
lat2=`gdalinfo geo_velocity.grd | sed -n '/^Upper Right/p' | awk '{print $5}' | sed 's/)//g'`
lim=-R$lon1/$lon2/$lat1/$lat2

if [[ $3 -eq -1 ]] && [[ $4 -eq 1 ]]; then
	cpt_min=`awk 'BEGIN {min = 65536} {if ($3+0 < min+0) min=$3} END {print min}' geo_velocity.xyz1`
	cpt_max=`awk 'BEGIN {max = -65536} {if ($3+0 > max+0) max=$3} END {print max}' geo_velocity.xyz1`
else
	cpt_min=$3
	cpt_max=$4
fi

gmt gmtset FORMAT_GEO_MAP ddd:mm:ssF
gmt gmtset MAP_ANNOT_OFFSET_PRIMARY 0.07i MAP_ANNOT_OBLIQUE 34 MAP_FRAME_WIDTH 0.08i MAP_SCALE_HEIGHT 0.1i MAP_LABEL_OFFSET 0.08i MAP_TICK_LENGTH_PRIMARY 0.1i PS_MEDIA B3 MAP_FRAME_TYPE fancy COLOR_NAN white FONT 10p

gmt makecpt -Cbgyr.cpt -T$cpt_min/$cpt_max/0.01 -I -V -Z > g.cpt
gmt psbasemap $lim -JM15 -Bxaf -Byaf -BWeSn -P -K -V  > geo_velocity.ps
gmt psxy geo_velocity.xyz1 $lim -JM15 -Cg.cpt -V -Sc0.01i -K -O >> geo_velocity.ps
gmt psscale -Cg.cpt -Dx4.1i/0.35i+w1.65i/0.125i+h -I -O -Baf >> geo_velocity.ps
gmt psconvert geo_velocity.ps -Tg -A -E800 -W+g
mv *.tiff $output_dir && rm *


cd $dir/geo && rm -rf PLOT2
mkdir -vp PLOT2 && cd PLOT2
mask.py ../geo_timeseries*.h5 -m ../geo_maskTempCoh.h5 -o ./geo_timeseries.h5
info.py geo_timeseries.h5 --date > date.tab
ref_date=`sed -n '1p' date.tab`
reference_date.py geo_timeseries.h5 --ref-date $ref_date
cat date.tab | while read a;do save_gmt.py geo_timeseries.h5 $a -o $a.disp;done

if [[ $5 -eq -1 ]] && [[ $6 -eq 1 ]]; then
	cpt_min_ts=$(echo "`info.py geo_timeseries.h5 --dset timeseries | sed -n '/^dataset min/p' | awk '{print $5}'` * 1000" | bc)
	cpt_max_ts=$(echo "`info.py geo_timeseries.h5 --dset timeseries | sed -n '/^dataset min/p' | awk '{print $7}'` * 1000" | bc)
else
	cpt_min_ts=$5
	cpt_max_ts=$6
fi

gmt gmtset FORMAT_GEO_MAP ddd:mm:ssF
gmt gmtset MAP_ANNOT_OFFSET_PRIMARY 0.07i MAP_ANNOT_OBLIQUE 34 MAP_FRAME_WIDTH 0.08i MAP_SCALE_HEIGHT 0.1i MAP_LABEL_OFFSET 0.08i MAP_TICK_LENGTH_PRIMARY 0.1i PS_MEDIA B3 MAP_FRAME_TYPE fancy COLOR_NAN white FONT 10p

gmt makecpt -Cbgyr.cpt -T$cpt_min_ts/$cpt_max_ts/0.01 -I -V -Z > g.cpt

ls *.disp > ts.tab
cat ts.tab | while read a;do
	gmt grdmath $a 0 NAN 1000 MUL = $a.grd
	gmt grd2xyz $a.grd > $a.xyz
	sed -e '/NaN/d' $a.xyz > $a.xyz1
	gmt psbasemap $lim -JM15 -Bxaf -Byaf -BWeSn -P -K -V  > $a.ps
	gmt psxy $a.xyz1 $lim -JM15 -Cg.cpt -V -Sc0.01i -K -O >> $a.ps
	gmt psscale -Cg.cpt -Dx4.1i/0.35i+w1.65i/0.125i+h -I -O -Baf >> $a.ps
	gmt psconvert $a.ps -Tg -A -E800 -W+g
	rm -f $a.ps $a.xyz $a.xyz1 $a.pgw $a.grd
done
ts0=`sed -n '1p' ts.tab`
gmt grdmath $ts0 1000 MUL = $ts0.grd
gmt grd2xyz $ts0.grd > $ts0.xyz
sed -e '/NaN/d' $ts0.xyz > $ts0.xyz1
gmt psbasemap $lim -JM15 -Bxaf -Byaf -BWeSn -P -K -V  > $ts0.ps
gmt psxy $ts0.xyz1 $lim -JM15 -Cg.cpt -V -Sc0.01i -K -O >> $ts0.ps
gmt psscale -Cg.cpt -Dx4.1i/0.35i+w1.65i/0.125i+h -I -O -Baf >> $ts0.ps
gmt psconvert $ts0.ps -Tg -A -E800 -W+g
mv *.tiff $output_dir && rm *

cd $output_dir
ls *.tiff > tiff_tab
cat tiff_tab | while read a;do gdaladdo -ro $a;done
rm tiff_tab
