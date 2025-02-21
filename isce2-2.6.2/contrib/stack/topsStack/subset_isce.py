#!/usr/bin/env python3

from mintpy.utils import utils as ut
from mintpy.utils import readfile, writefile, isce_utils
import argparse
import os
import shutil
import glob

class customArgparseFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    '''
    For better help message that also shows the defaults.
    '''
    pass

def cmdLineParse():
    '''
    Command Line Parser.
    '''
    parser = argparse.ArgumentParser(description='Make a subset from isce merged products.',
            formatter_class=customArgparseFormatter,
            epilog = '''

Example:

subset_isce.py -d /tmp/data/Xian/process -b '33.6 33.9 113.0 113.4' -r 10 -a 2

''')
    parser.add_argument('-d','--dir', type=str, required=True, help='work directory', dest='workdir')
    parser.add_argument('-b','--bbox',type=str, required=True, help="Lat/Lon Bounding SNWE. -- Example : '33.6 33.8 113.0 113.3'", dest='bbox')
    parser.add_argument('-r', '--range', type=int, required=True, default=1, help='Number of range looks. Default: 1', dest='rglooks')
    parser.add_argument('-a', '--azimuth', type=int, required=True, default=1, help='Number of azimuth looks. Default: 1', dest='azlooks')
    values = parser.parse_args()

    return values


def main(inps):
    '''
    The main driver.
    '''
    
    # create .rsc file for geom_files
    meta_files = sorted(glob.glob(os.path.join(inps.workdir, 'reference/IW*.xml')))
    if len(meta_files) > 0:
        meta_file = meta_files[0]
    metadata = {}
    rsc_file = os.path.join(os.path.dirname(meta_file), 'data.rsc')
    if os.path.isdir(os.path.join(inps.workdir, 'merged/geom_reference_full')):
        geom_dir = os.path.join(inps.workdir, 'merged/geom_reference_full')
    else:
        geom_dir = os.path.join(inps.workdir, 'merged/geom_reference')
    metadata = isce_utils.extract_isce_metadata(
            meta_file,
            geom_dir=geom_dir,
            rsc_file=rsc_file,
            update_mode=True)[0]
    
    geometry_prefixs = ['hgt', 'lat', 'lon', 'los', 'shadowMask', 'waterMask', 'incLocal']
    geom_files = [f'{i}.rdr' for i in geometry_prefixs]
    geom_files = [os.path.join(geom_dir, i) for i in geom_files]
    geom_files = [i for i in geom_files if os.path.isfile(i)]
    for geom_file in geom_files:
        # prepare metadata for current file
        geom_meta = {**metadata}
        if os.path.isfile(geom_file+'.xml'):
            geom_meta.update(readfile.read_attribute(geom_file, metafile_ext='.xml'))
        else:
            geom_meta.update(readfile.read_attribute(geom_file))
        # write .rsc file
        rsc_file = geom_file+'.rsc'
        writefile.write_roipac_rsc(geom_meta, rsc_file,
                                   update_mode=True,
                                   print_msg=True)
    
    # translate geo_box to sarpix
    if inps.bbox is not None:
        box = [float(val) for val in inps.bbox.split()]

    geo_box = (box[2], box[1], box[3], box[0])     # geo_box - tuple of 4 float, indicating the UL/LR lon/lat
    pix_box = None 
    lat = os.path.join(geom_dir, 'lat.rdr')
    lon = os.path.join(geom_dir, 'lon.rdr')
    dem = os.path.join(geom_dir, 'hgt.rdr')
    lookup_file = [lat, lon]
    atr = readfile.read_attribute(dem)
    coord = ut.coordinate(atr, lookup_file=lookup_file)
    pix_box = coord.bbox_geo2radar(geo_box)
    pix_box = coord.check_box_within_data_coverage(pix_box)
    print(f'input bounding box of interest in lalo: {geo_box}')
    print(f'box to subset for datasets in y/x: {pix_box}')
    
    # copy origin data to {dir}_full
    geom_path1 = os.path.join(inps.workdir, 'merged/geom_reference')
    geom_path2 = os.path.join(inps.workdir, 'merged/geom_reference_full')
    if not os.path.isdir(geom_path2):
        shutil.copytree(geom_path1, geom_path2)
    ifg_path1 = os.path.join(inps.workdir, 'merged/interferograms')
    ifg_path2 = os.path.join(inps.workdir, 'merged/interferograms_full')
    if not os.path.isdir(ifg_path2):
        shutil.copytree(ifg_path1, ifg_path2)
    
    # subset geom datasets from {dir}_full to {dir}
    geo_files = [f'{i}.rdr' for i in geometry_prefixs]
    geo_files = [os.path.join(geom_path1, i) for i in geo_files]
    geo_files = [i for i in geo_files if os.path.isfile(i)]
    geo_files_full = [f'{i}.rdr' for i in geometry_prefixs]
    geo_files_full = [os.path.join(geom_path2, i) for i in geo_files_full]
    geo_files_full = [i for i in geo_files_full if os.path.isfile(i)]
    shutil.rmtree(geom_path1)
    for i in range(len(geo_files_full)):
        os.system('subset.py %s -y %d %d -x %d %d -o %s' % (geo_files_full[i], pix_box[1], pix_box[3], pix_box[0], pix_box[2], geo_files[i]))
    
    # subset ifg datasets from {dir}_full to {dir}
    ifg_files = [os.path.join(inps.workdir, 'merged/interferograms', i, 'fine.int') for i in os.listdir(os.path.join(inps.workdir, 'merged/interferograms_full'))]
    ifg_files_full = [os.path.join(inps.workdir, 'merged/interferograms_full', i, 'fine.int') for i in os.listdir(os.path.join(inps.workdir, 'merged/interferograms_full'))]
    ifg_files_full = [i for i in ifg_files_full if os.path.isfile(i)]
    shutil.rmtree(os.path.join(inps.workdir, 'merged/interferograms'))
    for i in range(len(ifg_files_full)):
        os.system('subset.py %s -y %d %d -x %d %d -o %s' % (ifg_files_full[i], pix_box[1], pix_box[3], pix_box[0], pix_box[2], ifg_files[i]))
        
    # calculate ave.mli and subset in {dir}/merged/SLC
    os.chdir(os.path.join(inps.workdir, 'merged/SLC'))
    os.system("ls -d */ > tab && sed -i 's/\///g' tab")
    os.system('cat tab | while read a;do multi_look.py -i $a/$a.slc.full.vrt -o $a/$a.mli -r %d -a %d;done' % (inps.rglooks, inps.azlooks))
    os.system('''find */ -name "*.mli" > mli.list''')
    os.system('image_ave.py -i mli.list -o ave_full.mli')
    os.system('subset.py ave_full.mli -y %d %d -x %d %d -o ave.mli' % (pix_box[1], pix_box[3], pix_box[0], pix_box[2]))
    os.system('rm tab mli.list ave_full.mli*')
    
    return 0


if __name__ == '__main__':
    '''
    Makes the script executable.
    '''

    inps = cmdLineParse()
    main(inps)
