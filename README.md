InSAR_process
Open Source SBAS-InSAR Framework From Comma "douhao_zy@163.com"

构建镜像前，请先填充config_files/cdsapirc, config_files/netrc中的账户信息！
docker build --rm --force-rm -t insar_process:v1.0 -f Dockerfile .
xhost +
docker run -itd --privileged=true --net=host --shm-size 4g --gpus=all --name insar -e DISPLAY=$DISPLAY -v /media/comma/InSAR2/Mexico/:/opt/data insar_process:v1.0
docker exec -it insar /bin/bash

进入容器...
runfiles_path=/opt/data/descending/process/run_files
dem_directory=/opt/data/descending/DEM
nasadem_directory=/opt/data/nasaDEM
dem_bbox="17 20 -101 -98"
working_directory=/opt/data/descending/process
slc_directory=/opt/data/descending/SLC
aux_directory=/opt/data/descending/AuxDir
orbit_directory=/opt/data/Orbits
slc_bbox='19.3 19.5 -99.15 -98.86'
azimuth_looks=2
range_looks=10
reference_date=20201206
num_connections=3
filter_strength=0.4
nfft=32
swath_num='3'
cc_thres=0.4
pwr_thres=0.1
num_process=2
num_process4topo=1
snr_misreg_threshold=10
num_overlap_connections=3
esd_coherence_threshold=0.85
parallel_num=4
subset=1

mkdir -vp $dem_directory && cd $dem_directory && dem.py -a stitch -b $dem_bbox -r -t nasadem -c -l -d $nasadem_directory && fixImageXml.py -i *.dem.wgs84 -f && rm -f $nasadem_directory/demLat* $nasadem_directory/*.hgt

mkdir -vp $working_directory $slc_directory $aux_directory $orbit_directory && dem=`find $dem_directory -name "*.dem.wgs84"` && stackSentinel.py --working_directory="$working_directory" --slc_directory="$slc_directory" --dem="$dem" --aux_directory="$aux_directory" --orbit_directory="$orbit_directory" --bbox="$slc_bbox" --azimuth_looks="$azimuth_looks" --range_looks="$range_looks" --reference_date="$reference_date" --num_connections="$num_connections" --filter_strength="$filter_strength" --nfft="$nfft" --swath_num="$swath_num" --cc_thres="$cc_thres" --pwr_thres="$pwr_thres" --num_process="$num_process" --num_process4topo="$num_process4topo" --snr_misreg_threshold="$snr_misreg_threshold" --num_overlap_connections="$num_overlap_connections" --esd_coherence_threshold="$esd_coherence_threshold" --useGPU

cd $runfiles_path && split_file.sh $parallel_num
chmod 777 -R $runfiles_path/*

${runfiles_path}/run_01*
for ((idx=1; idx<=$(($(ls ${runfiles_path}/run_02* 2>/dev/null | wc -l)-1)); idx++)); do ${runfiles_path}/run_02*_${idx}; done
${runfiles_path}/run_03*
${runfiles_path}/run_04*
for ((idx=1; idx<=$(($(ls ${runfiles_path}/run_05* 2>/dev/null | wc -l)-1)); idx++)); do ${runfiles_path}/run_05*_${idx}; done
for ((idx=1; idx<=$(($(ls ${runfiles_path}/run_06* 2>/dev/null | wc -l)-1)); idx++)); do ${runfiles_path}/run_06*_${idx}; done
for ((idx=1; idx<=$(($(ls ${runfiles_path}/run_07* 2>/dev/null | wc -l)-1)); idx++)); do ${runfiles_path}/run_07*_${idx}; done
${runfiles_path}/run_08*
for ((idx=1; idx<=$(($(ls ${runfiles_path}/run_09* 2>/dev/null | wc -l)-1)); idx++)); do ${runfiles_path}/run_09*_${idx}; done
for ((idx=1; idx<=$(($(ls ${runfiles_path}/run_10* 2>/dev/null | wc -l)-1)); idx++)); do ${runfiles_path}/run_10*_${idx}; done
${runfiles_path}/run_11*
${runfiles_path}/run_12*
for ((idx=1; idx<=$(($(ls ${runfiles_path}/run_13* 2>/dev/null | wc -l)-1)); idx++)); do ${runfiles_path}/run_13*_${idx}; done
for ((idx=1; idx<=$(($(ls ${runfiles_path}/run_14* 2>/dev/null | wc -l)-1)); idx++)); do ${runfiles_path}/run_14*_${idx}; done
if [ "${subset}" -eq 1 ]; then subset_isce.py -d ${working_directory} -b "$slc_bbox" -r ${range_looks} -a ${azimuth_looks}; fi
for ((idx=1; idx<=$(($(ls ${runfiles_path}/run_15* 2>/dev/null | wc -l)-1)); idx++)); do ${runfiles_path}/run_15*_${idx}; done
pre_unwrap.sh ${working_directory} ${range_looks} ${azimuth_looks}
for ((idx=1; idx<=$(($(ls ${runfiles_path}/run_16* 2>/dev/null | wc -l)-1)); idx++)); do ${runfiles_path}/run_16*_${idx}; done
unw=`find ${working_directory} -name "filt_fine.unw" | sed -n "1p"` && gdal_translate ${unw} -b 1 -of ISCE ${working_directory}/magnitude && generate_mask.py ${working_directory}/magnitude --min 1 -o ${working_directory}/waterMask.h5 && rm ${working_directory}/magnitude*

smallbaselineApp.py ./mexicoSentinel.txt
退出容器...

docker stop insar && docker rm insar
