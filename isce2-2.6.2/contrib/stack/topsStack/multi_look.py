#!/usr/bin/env python3

from osgeo import gdal
import numpy as np
import warnings
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
    parser = argparse.ArgumentParser(description='Take integer number of looks.',
            formatter_class=customArgparseFormatter,
            epilog = '''

Example:

multi_look.py -i input.file -o output.file -r 10 -a 2

''')
    parser.add_argument('-i','--input', type=str, required=True, help='Input ISCEproduct with a corresponding .xml file.', dest='infile')
    parser.add_argument('-o','--output',type=str, default=None, help='Output ISCE DEproduct with a corresponding .xml file.', dest='outfile')
    parser.add_argument('-r', '--range', type=int, default=1, help='Number of range looks. Default: 1', dest='rglooks')
    parser.add_argument('-a', '--azimuth', type=int, default=1, help='Number of azimuth looks. Default: 1', dest='azlooks')

    values = parser.parse_args()
    if (values.rglooks == 1) and (values.azlooks == 1):
        print('Nothing to do. One look requested in each direction. Exiting ...')
        sys.exit(0)

    return values


def multi_look(infile, outfile, nlook_r, nlook_c, n_valid_thre = 0.5):
	ds = gdal.Open(infile, gdal.GA_ReadOnly)
	data = ds.GetRasterBand(1).ReadAsArray()
	intensity = np.square(data.real) + np.square(data.imag)
	length, width = intensity.shape
	length_ml = int(np.floor(length/nlook_c))
	width_ml = int(np.floor(width/nlook_r))
	intensity_reshape = intensity[:length_ml*nlook_c,:width_ml*nlook_r].reshape(length_ml, nlook_c, width_ml, nlook_r)
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', RuntimeWarning)
		intensity_ml = np.nanmean(intensity_reshape, axis=(1, 3))
	n_valid = np.sum(~np.isnan(intensity_reshape), axis=(1, 3))
	bool_invalid = n_valid < n_valid_thre*nlook_r*nlook_c
	intensity_ml[bool_invalid] = np.nan

	driver = gdal.GetDriverByName('ISCE')
	mli=driver.Create(outfile, width_ml, length_ml , 1, gdal.GDT_Float32)
	mli.GetRasterBand(1).WriteArray(intensity_ml)
	del mli


def main(inps):
    '''
    The main driver.
    '''
    
    if inps.infile.endswith('.xml'):
        inFileXml = inps.infile
        inFile = os.path.splitext(inps.infile)[0]
    else:
        inFile = inps.infile
        inFileXml = inps.infile + '.xml'
        
    if inps.outfile is None:
        spl = os.path.splitext(inFile)
        ext = '.{0}alks_{1}rlks'.format(inps.azlooks, inps.rglooks)
        outFile = spl[0] + ext + spl[1]
    elif inps.outfile.endswith('.xml'):
        outFile = os.path.splitext(inps.outfile)[0]
    else:
        outFile = inps.outfile
    
    print('Output filename : {0}'.format(outFile))
    multi_look(inFile, outFile, inps.rglooks, inps.azlooks)
    
    return outFile
        

if __name__ == '__main__':
    '''
    Makes the script executable.
    '''

    inps = cmdLineParse()
    main(inps)
