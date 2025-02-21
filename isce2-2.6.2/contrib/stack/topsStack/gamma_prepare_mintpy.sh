reference=20201206
rlks=10

mkdir geometry interferograms mintpy
cp outputdem/"$reference".lookup_tabke_final geometry/sim_"$reference"_"$rlks"rlks.UTM_TO_RDC
cp outputdem/"$reference".diff.geo.par geometry/sim_"$reference"_"$rlks"rlks.diff_par
cp outputdem/"$reference".dem.final geometry/sim_"$reference"_"$rlks"rlks.rdc.dem
cp outputdem/"$reference".dem.cropped.par geometry/sim_"$reference"_"$rlks"rlks.utm.dem.par

cd interferograms/
mkdir `cat ../bperp.ascii.all | awk '{print $2"_"$3}'`
cd ..
cat bperp.ascii.all | awk '{print $2"_"$3}' > interferograms.tab

cat interferograms.tab | while read a;do
	cp diff2d/$a.adf.cc interferograms/$a/filt_"$a"_"$rlks"rlks.cor
	cp diff2d/$a.adf.unw interferograms/$a/diff_filt_"$a"_"$rlks"rlks.unw
	cp diff2d/$a.base interferograms/$a/"$a"_"$rlks"rlks.baseline
	cp diff2d/$a.off interferograms/$a/"$a"_"$rlks"rlks.off
done

run_all.pl bperp.ascii.all 'cp mli/$2.rmli.par interferograms/$2_$3/$2_10rlks.ramp.par'
run_all.pl bperp.ascii.all 'cp mli/$3.rmli.par interferograms/$2_$3/$3_10rlks.ramp.par'
run_all.pl bperp.ascii.all 'SLC_corners mli/$2.rmli.par - - > interferograms/$2_$3/$2_10rlks.ramp.corner'
run_all.pl bperp.ascii.all 'sed -n -i '3,6p' interferograms/$2_$3/$2_10rlks.ramp.corner'
cat interferograms.tab | while read a;do awk '{print $3,$6}' interferograms/$a/*_10rlks.ramp.corner > interferograms/$a/tmp;done
run_all.pl bperp.ascii.all 'mv interferograms/$2_$3/tmp interferograms/$2_$3/$2_10rlks.ramp.corner'
run_all.pl bperp.ascii.all 'base_perp diff2d/$2_$3.base RSLC/$2.rslc.par diff2d/$2_$3.off > interferograms/$2_$3/$2_$3_10rlks.base_perp'

cd mintpy
smallbaselineApp.py gamma.template
