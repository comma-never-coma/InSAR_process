#! /bin/bash
##########################################################
# Script that:
# -generate ave.mli and ave.cor for unwrapping
#       
#########################################################################################################

##########################################################################################################
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo " `basename $0`: generate ave.mli and ave.cor for unwrapping "
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
if [ $# -ne 3 ]; then
        echo " pre_unwrap.sh <work_directory> [rlooks] [alooks]"
        echo "      work_directory: ISCE process path. "
        echo "      rlooks: Number of range looks. "
        echo "      alooks: Number of azimuth looks. "
        exit 1
fi

work_dir=$1
rlooks=$2
alooks=$3

if [ ! -f "$1/merged/SLC/ave.mli" ];then
    cd $1/merged/SLC 
    ls -d */ > tab && sed -i 's/\///g' tab
    cat tab | while read a;do multi_look.py -i $a/$a.slc.full.vrt -o $a/$a.mli -r $2 -a $3;done 
    find */ -name "*.mli" > mli.list
    image_ave.py -i mli.list -o ave.mli 
    rm tab mli.list
fi

cd $1/merged/interferograms
find */ -name "filt_fine.cor" > cor.list
image_ave.py -i cor.list -o ave.cor
rm cor.list
