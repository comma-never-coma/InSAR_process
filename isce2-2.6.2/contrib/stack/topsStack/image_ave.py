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
    parser = argparse.ArgumentParser(description='Take average of images.',
            formatter_class=customArgparseFormatter,
            epilog = '''

Example:

image_ave.py -i inputfile.list -o outputfile

''')
    parser.add_argument('-i','--input', type=str, required=True, help='Inputfile list consist of ISCEproduct with a corresponding .xml file.', dest='infile')
    parser.add_argument('-o','--output',type=str, required=True, help='Output ISCE DEproduct with a corresponding .xml file.', dest='outfile')

    values = parser.parse_args()

    return values


def image_ave(infile, outfile):
	with open(infile, 'r') as f:
		filelist = f.read()
	files = filelist.split('\n')
	width = gdal.Open(files[0], gdal.GA_ReadOnly).GetRasterBand(1).XSize
	length = gdal.Open(files[0], gdal.GA_ReadOnly).GetRasterBand(1).YSize
	mli = np.zeros((len(files)-1, length, width))
	for i in range(len(files)-1):
		ds = gdal.Open(os.path.abspath(files[i]), gdal.GA_ReadOnly)
		data = ds.GetRasterBand(1).ReadAsArray()
		mli[i] = data
	mli_ave = np.nanmean(mli, axis=0)
	
	driver = gdal.GetDriverByName('ISCE')
	ave_mli = driver.Create(outfile, width, length , 1, gdal.GDT_Float32)
	ave_mli.GetRasterBand(1).WriteArray(mli_ave)
	del ave_mli
		

def main(inps):
    '''
    The main driver.
    '''
    
    print('Output filename : {0}'.format(inps.outfile))
    image_ave(inps.infile, inps.outfile)
    
    return inps.outfile


if __name__ == '__main__':
    '''
    Makes the script executable.
    '''

    inps = cmdLineParse()
    main(inps)
