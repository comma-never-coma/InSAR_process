#! /bin/bash
##########################################################
# Script that:
# -combine asc and dsc mintpy results to 2D deformation (EW & UD)
#       
#########################################################################################################

##########################################################################################################
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo " `basename $0`: run_msbas.sh  "
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
if [ $# -ne 5 ]; then
        echo " run_msbas.sh <asc_mintpy_dir> <dsc_mintpy_dir> [bbox] [step] <output_dir>"
        echo "      asc_mintpy_dir: ascending data mintpy path "
        echo "      dsc_mintpy_dir: descending data mintpy path "
        echo "      bbox: area of interest, eg: '37.533027778 39.325386 108.404008333 110.929355556' "
        echo "      step: latitude and longitude step, eg: 0.000925926(100m), 0.000833334(90m), 0.000555556(60m), 0.000462963(50m), 0.000277778(30m), 0.000185185(20m), 0.000092593(10m)"
        echo "      output_dir: output data process path "
        exit 1
fi


asc_mintpy_dir=$1
dsc_mintpy_dir=$2
bbox=$3
step=$4
output_dir=$5


cd $asc_mintpy_dir/geo
sed -n '6,$p' ../coherenceSpatialAvg.txt | awk '{print $1}' > int.tab
asc_timeseries=($(ls geo_timeseries*.h5))
# 初始化最新文件信息
latest_asc=""
latest_asc_timestamp=0
# 遍历数组，比较文件时间戳以找到最新文件
for file in "${asc_timeseries[@]}"; do
  current_asc_timestamp=$(stat -c %Y "$file")
  if [ "$current_asc_timestamp" -gt "$latest_asc_timestamp" ]; then
    latest_asc_timestamp=$current_asc_timestamp
    latest_asc=$file
  fi
done
mask.py $latest_asc -m geo_maskTempCoh.h5 -o geo_timeseries_asc.h5
cat int.tab | while read a;do save_gdal.py geo_timeseries_asc.h5 -d $a;done
sed -i 's/_/ /g' int.tab
awk 'NR>5 {print $4}' ../coherenceSpatialAvg.txt > ./Bperp.txt
mkdir -p asc
mv *.tif ./asc
ls asc/*.tif > asc_tif.tab
paste asc_tif.tab Bperp.txt int.tab > asc.txt
rm asc_tif.tab Bperp.txt int.tab
mv asc.txt asc geo_timeseries_asc.h5 $output_dir


cd $dsc_mintpy_dir/geo
sed -n '6,$p' ../coherenceSpatialAvg.txt | awk '{print $1}' > int.tab
dsc_timeseries=($(ls geo_timeseries*.h5))
# 初始化最新文件信息
latest_dsc=""
latest_dsc_timestamp=0
# 遍历数组，比较文件时间戳以找到最新文件
for file in "${dsc_timeseries[@]}"; do
  current_dsc_timestamp=$(stat -c %Y "$file")
  if [ "$current_dsc_timestamp" -gt "$latest_dsc_timestamp" ]; then
    latest_dsc_timestamp=$current_dsc_timestamp
    latest_dsc=$file
  fi
done
mask.py $latest_dsc -m geo_maskTempCoh.h5 -o geo_timeseries_dsc.h5
cat int.tab | while read a;do save_gdal.py geo_timeseries_dsc.h5 -d $a;done
sed -i 's/_/ /g' int.tab
awk 'NR>5 {print $4}' ../coherenceSpatialAvg.txt > ./Bperp.txt
mkdir -p dsc
mv *.tif ./dsc
ls dsc/*.tif > dsc_tif.tab
paste dsc_tif.tab Bperp.txt int.tab > dsc.txt
rm dsc_tif.tab Bperp.txt int.tab
mv dsc.txt dsc geo_timeseries_dsc.h5 $output_dir


cd $output_dir
#############################mintpy: geocode#############################
ls asc/*.tif > asc.tab
ls dsc/*.tif > dsc.tab
sed -i 's/.tif//g' asc.tab
sed -i 's/.tif//g' dsc.tab
aoi=($bbox)
cat asc.tab | while read a;do gdal_translate -projwin ${aoi[2]} ${aoi[1]} ${aoi[3]} ${aoi[0]} -tr $step -$step -of GTiff $a.tif $a.crop.tif;done
cat dsc.tab | while read a;do gdal_translate -projwin ${aoi[2]} ${aoi[1]} ${aoi[3]} ${aoi[0]} -tr $step -$step -of GTiff $a.tif $a.crop.tif;done
rm asc.tab dsc.tab
sed -i 's/.tif/.crop.tif/g' asc.txt
sed -i 's/.tif/.crop.tif/g' dsc.txt
#########################################################################
tif1=`cat asc.txt | awk 'NR==1 {print$1}'`
width1=`info.py $tif1 | grep WIDTH | awk '{print $2}'`
length1=`info.py $tif1 | grep LENGTH | awk '{print $2}'`
tif2=`cat dsc.txt | awk 'NR==1 {print$1}'`
width2=`info.py $tif2 | grep WIDTH | awk '{print $2}'`
length2=`info.py $tif2 | grep LENGTH | awk '{print $2}'`

if [ $width1 -lt $width2 ]
then
	width=$width1
else
	width=$width2
fi

if [ $length1 -lt $length2 ]
then
	length=$length1
else
	length=$length2
fi

inc1=`info.py geo_timeseries_asc.h5 | grep CENTER_INCIDENCE_ANGLE | awk '{print $2}'`
azi1=`info.py geo_timeseries_asc.h5 | grep HEADING | awk '{print $2}'`
inc2=`info.py geo_timeseries_dsc.h5 | grep CENTER_INCIDENCE_ANGLE | awk '{print $2}'`
azi2=`info.py geo_timeseries_dsc.h5 | grep HEADING | awk '{print $2}'`
echo "FORMAT=2" > header.txt
echo "FILE_SIZE=$width,$length" >> header.txt
echo "C_FLAG=0" >> header.txt
echo "R_FLAG=0" >> header.txt
echo "I_FLAG=0" >> header.txt
echo "SET=0,000000,$azi1,$inc1,asc.txt" >> header.txt
echo "SET=0,000000,$azi2,$inc2,dsc.txt" >> header.txt
msbas header.txt
