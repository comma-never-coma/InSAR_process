#!/usr/bin/env python3

from osgeo import gdal
import numpy as np
import argparse
import os

class customArgparseFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    '''
    For better help message that also shows the defaults.
    '''
    pass

def cmdLineParse():
    '''
    Command Line Parser.
    '''
    parser = argparse.ArgumentParser(description='Select pslike points based on coherence.',
            formatter_class=customArgparseFormatter,
            epilog = '''

Example:

select_pslike.py -t 0.5 -p 0.8 -r 16 -n 4

''')
    parser.add_argument('-t','--threshold', type=float, required=True, help='Coherence threshold to be selected(0~1).', dest='threshold')
    parser.add_argument('-p','--percent',type=float, required=True, help='Coherenct coverage on the image during acquisitions.', dest='percent')
    parser.add_argument('-r','--ram', type=int, required=True, help='Max amount of memory in GB to use.', dest='ram')
    parser.add_argument('-n','--numprocess',type=int, required=True, help='Number of processors to be used for calculation.', dest='numprocess')

    values = parser.parse_args()

    return values

def select_pslike(threshold, percent, ram, numprocess):
	os.system('find . -name "fine.cor" > cor.list')
	os.system("sed -i 's/\/fine.cor//g' cor.list")
	os.system('cat cor.list | while read a;do geocode.py $a/fine.cor -d band2 --lat-file ../geom_reference/lat.rdr --lon-file ../geom_reference/lon.rdr --outdir $a --ram %d -n %d;done' % (ram, numprocess))
	os.system('find . -name "geo_band2.cor" > cor.geo.list')

	with open('cor.geo.list', 'r') as f:
		filelist = f.read()
	files = filelist.split('\n')
	width = gdal.Open(files[0]+'.vrt', gdal.GA_ReadOnly).GetRasterBand(1).XSize
	length = gdal.Open(files[0]+'.vrt', gdal.GA_ReadOnly).GetRasterBand(1).YSize
	cc = np.zeros((length, width, len(files)-1))
	for i in range(len(files)-1):
		ds = gdal.Open(os.path.abspath(files[i]+'.vrt'), gdal.GA_ReadOnly)
		data = ds.GetRasterBand(1).ReadAsArray()
		cc[:,:,i] = data

	k = cc.shape[2]
	IDX = np.arange(length*width)
	coh = np.reshape(cc, [length*width, k])
	coh = coh >= threshold
	coh = np.sum(coh, axis=1)
	num = np.where(coh<k*percent)
	IDX = np.delete(IDX, num)
	coh = np.delete(coh, num)
	[Y,X] = np.unravel_index(IDX, [length, width])

	with open('X','w') as f:
		np.savetxt(f,X,fmt='%d',delimiter='\t')
	with open('Y','w') as f:
		np.savetxt(f,Y,fmt='%d',delimiter='\t')
	os.system('paste X Y > tcp_ps')
	os.system('rm X Y cor.list cor.geo.list')

	trans = gdal.Open(files[0], gdal.GA_ReadOnly).GetGeoTransform()
	px = trans[0] + X * trans[1] + Y * trans[2]
	py = trans[3] + X * trans[4] + Y * trans[5]
	with open('px','w') as f:
		np.savetxt(f,px,fmt='%.15f',delimiter='\t')
	with open('py','w') as f:
		np.savetxt(f,py,fmt='%.15f',delimiter='\t')
	os.system('paste px py > geo_tcp')
	os.system('rm px py')

def main(inps):
    '''
    The main driver.
    '''

    select_pslike(inps.threshold, inps.percent, inps.ram, inps.numprocess)
    
    return 0


if __name__ == '__main__':
    '''
    Makes the script executable.
    '''

    inps = cmdLineParse()
    main(inps)
