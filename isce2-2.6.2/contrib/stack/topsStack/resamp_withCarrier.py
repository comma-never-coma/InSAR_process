#!/usr/bin/env python3

import isce
import isceobj
import stdproc
from stdproc.stdproc import crossmul
import numpy as np
from isceobj.Util.Poly2D import Poly2D
import argparse
import os
import copy
import s1a_isce_utils as ut
from isceobj.Sensor.TOPS import createTOPSSwathSLCProduct


def createParser():
    parser = argparse.ArgumentParser( description='Resampling burst by burst SLCs ')

    parser.add_argument('-m', '--reference', dest='reference', type=str, required=True,
            help='Directory with reference acquisition')

    parser.add_argument('-s', '--secondary', dest='secondary', type=str, required=True,
            help='Directory with secondary acquisition')

    parser.add_argument('-o', '--coregdir', dest='coreg', type=str, default='coreg_secondary',
            help='Directory with coregistered SLCs and IFGs')

    parser.add_argument('-a', '--azimuth_misreg', dest='misreg_az', type=str, default=0.0,
            help='File name with the azimuth misregistration')

    parser.add_argument('-r', '--range_misreg', dest='misreg_rng', type=str, default=0.0,
            help='File name with the range misregistration')

    parser.add_argument('--noflat', dest='noflat', action='store_true', default=False,
            help='To turn off flattening. False: flattens the SLC. True: turns off flattening.')

    parser.add_argument('-v', '--overlap', dest='overlap', action='store_true', default=False,
            help='Is this an overlap burst slc. default: False')

    parser.add_argument('-d', '--overlapDir', dest='overlapDir', type=str, default='overlap',
            help='reference overlap directory')
            
    parser.add_argument('-useGPU', '--useGPU', dest='useGPU',action='store_true', default=False,
            help='Allow App to use GPU when available')

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)

def resampSecondaryCPU(mas, slv, rdict, outname, flatten):
    '''
    Resample burst by burst.
    '''

    azpoly = rdict['azpoly']
    rgpoly = rdict['rgpoly']
    azcarrpoly = rdict['carrPoly']
    dpoly = rdict['doppPoly']

    rngImg = isceobj.createImage()
    rngImg.load(rdict['rangeOff'] + '.xml')
    rngImg.setAccessMode('READ')

    aziImg = isceobj.createImage()
    aziImg.load(rdict['azimuthOff'] + '.xml')
    aziImg.setAccessMode('READ')

    inimg = isceobj.createSlcImage()
    inimg.load(slv.image.filename + '.xml')
    inimg.setAccessMode('READ')
    

    rObj = stdproc.createResamp_slc()
    rObj.slantRangePixelSpacing = slv.rangePixelSize
    rObj.radarWavelength = slv.radarWavelength
    rObj.azimuthCarrierPoly = azcarrpoly
    rObj.dopplerPoly = dpoly

    rObj.azimuthOffsetsPoly = azpoly
    rObj.rangeOffsetsPoly = rgpoly
    rObj.imageIn = inimg

    width = mas.numberOfSamples
    length = mas.numberOfLines
    imgOut = isceobj.createSlcImage()
    imgOut.setWidth(width)
    imgOut.filename = outname
    imgOut.setAccessMode('write')

    rObj.outputWidth = width
    rObj.outputLines = length
    rObj.residualRangeImage = rngImg
    rObj.residualAzimuthImage = aziImg
    rObj.flatten = flatten
    print(rObj.flatten)
    rObj.resamp_slc(imageOut=imgOut)

    imgOut.renderHdr()
    imgOut.renderVRT()
    return imgOut

def convertPoly2D(poly):
    '''
    Convert a isceobj.Util.Poly2D {poly} into zerodop.GPUresampslc.GPUresampslc.PyPloy2d
    '''
    from zerodop.GPUresampslc.GPUresampslc import PyPoly2d
    import itertools

    # get parameters from poly
    azimuthOrder = poly.getAzimuthOrder()
    rangeOrder = poly.getRangeOrder()
    azimuthMean = poly.getMeanAzimuth()
    rangeMean = poly.getMeanRange()
    azimuthNorm = poly.getNormAzimuth()
    rangeNorm = poly.getNormRange()

    # create the PyPoly2d object
    pPoly = PyPoly2d(azimuthOrder, rangeOrder, azimuthMean, rangeMean, azimuthNorm, rangeNorm)
    # copy the coeffs, need to flatten into 1d list
    pPoly.coeffs = list(itertools.chain.from_iterable(poly.getCoeffs()))

    # all done
    return pPoly

def resampSecondaryGPU(mas, slv, rdict, outname, flatten):
    '''
    Resample burst by burst with GPU.
    '''

    # import the GPU module
    from zerodop.GPUresampslc.GPUresampslc import PyResampSlc

    # get Poly2D objects from rdict and convert them into PyPoly2d objects
    azpoly = convertPoly2D(rdict['azpoly'])
    rgpoly = convertPoly2D(rdict['rgpoly'])
    azcarrpoly = convertPoly2D(rdict['carrPoly'])
    dpoly = convertPoly2D(rdict['doppPoly'])

    rngImg = isceobj.createImage()
    rngImg.load(rdict['rangeOff'] + '.xml')
    rngImg.setCaster('read', 'FLOAT')
    rngImg.setAccessMode('READ')
    rngImg.createImage()

    aziImg = isceobj.createImage()
    aziImg.load(rdict['azimuthOff'] + '.xml')
    aziImg.setCaster('read', 'FLOAT')
    aziImg.setAccessMode('READ')
    aziImg.createImage()
    
    inimg = isceobj.createSlcImage()
    inimg.load(slv.image.filename + '.xml')
    inimg.setAccessMode('READ')
    inimg.createImage()

    # create a GPU resample processor
    rObj = PyResampSlc()
    # set parameters
    rObj.slr = slv.rangePixelSize
    rObj.wvl = slv.radarWavelength
    # set polynomials
    rObj.azCarrier = azcarrpoly
    rObj.dopplerPoly = dpoly
    rObj.azOffsetsPoly = azpoly
    rObj.rgOffsetsPoly = rgpoly
    # need to create an empty rgCarrier poly
    rgCarrier = Poly2D()
    rgCarrier.initPoly(rangeOrder=0, azimuthOrder=0, coeffs=[[0.]])
    rgCarrier = convertPoly2D(rgCarrier)
    rObj.rgCarrier = rgCarrier
    
    # input secondary image
    rObj.slcInAccessor = inimg.getImagePointer()
    rObj.inWidth = inimg.getWidth()
    rObj.inLength = inimg.getLength()
    
    ####Setting reference values
    rObj.r0 = slv.startingRange
    rObj.refr0 = mas.startingRange
    rObj.refslr = mas.rangePixelSize
    rObj.refwvl = mas.radarWavelength
    
    # set output image
    width = mas.numberOfSamples
    length = mas.numberOfLines
    imgOut = isceobj.createSlcImage()
    imgOut.setWidth(width)
    imgOut.filename = outname
    imgOut.setAccessMode('write')
    imgOut.createImage()
    rObj.slcOutAccessor = imgOut.getImagePointer()
    
    rObj.outWidth = width
    rObj.outLength = length
    rObj.residRgAccessor = rngImg.getImagePointer()
    rObj.residAzAccessor = aziImg.getImagePointer()
    
    rObj.flatten = flatten
    print(rObj.flatten)
    
    # need to specify data type, only complex is currently supported
    rObj.isComplex = (inimg.dataType == 'CFLOAT')
    # run resampling
    rObj.resamp_slc()

    # finalize images
    inimg.finalizeImage()
    imgOut.finalizeImage()
    rngImg.finalizeImage()
    aziImg.finalizeImage()

    imgOut.renderHdr()
    imgOut.renderVRT()
    return imgOut

def main(iargs=None):
    '''
    Create coregistered overlap secondarys.
    '''
    inps = cmdLineParse(iargs)
    
    # see if the user compiled isce with GPU enabled
    run_GPU = False
    try:
        from zerodop.GPUresampslc.GPUresampslc import PyResampSlc
        run_GPU = True
    except:
        pass
        
    if inps.useGPU and not run_GPU:
        print("GPU mode requested but no GPU ISCE code found")
    
    # setting the respective version of geo2rdr for CPU and GPU
    if run_GPU and inps.useGPU:
        print('Using GPU for resampslc')
        resampSecondary = resampSecondaryGPU
    else:
        resampSecondary = resampSecondaryCPU
    
    referenceSwathList = ut.getSwathList(inps.reference)
    secondarySwathList = ut.getSwathList(inps.secondary)

    swathList = list(sorted(set(referenceSwathList+secondarySwathList)))

    for swath in swathList:
    
        ####Load secondary metadata
        reference = ut.loadProduct( os.path.join(inps.reference , 'IW{0}.xml'.format(swath)))
        secondary = ut.loadProduct( os.path.join(inps.secondary , 'IW{0}.xml'.format(swath)))
        if inps.overlap:
            referenceTop = ut.loadProduct(os.path.join(inps.reference, inps.overlapDir , 'IW{0}_top.xml'.format(swath)))
            referenceBottom = ut.loadProduct(os.path.join(inps.reference, inps.overlapDir , 'IW{0}_bottom.xml'.format(swath)))

    
        dt = secondary.bursts[0].azimuthTimeInterval
        dr = secondary.bursts[0].rangePixelSize

        if os.path.exists(str(inps.misreg_az)):
             with open(inps.misreg_az, 'r') as f:
                misreg_az = float(f.readline())
        else:
             misreg_az = 0.0

        if os.path.exists(str(inps.misreg_rng)):
             with open(inps.misreg_rng, 'r') as f:
                misreg_rg = float(f.readline())
        else:
             misreg_rg = 0.0

        ###Output directory for coregistered SLCs
        if not inps.overlap:
            outdir = os.path.join(inps.coreg,'IW{0}'.format(swath))
            offdir = os.path.join(inps.coreg,'IW{0}'.format(swath)) 
        else:
            outdir = os.path.join(inps.coreg, inps.overlapDir, 'IW{0}'.format(swath))
            offdir = os.path.join(inps.coreg, inps.overlapDir, 'IW{0}'.format(swath))
        os.makedirs(outdir, exist_ok=True)

    
        ####Indices w.r.t reference
        burstoffset, minBurst, maxBurst = reference.getCommonBurstLimits(secondary)
        secondaryBurstStart = minBurst +  burstoffset
        secondaryBurstEnd = maxBurst
    
        relShifts = ut.getRelativeShifts(reference, secondary, minBurst, maxBurst, secondaryBurstStart)
        if inps.overlap:
            maxBurst = maxBurst - 1 ###For overlaps
    
        print('Shifts: ', relShifts)
    
    ####Can corporate known misregistration here
    
        apoly = Poly2D()
        apoly.initPoly(rangeOrder=0,azimuthOrder=0,coeffs=[[0.]])
    
        rpoly = Poly2D()
        rpoly.initPoly(rangeOrder=0,azimuthOrder=0,coeffs=[[0.]])

    
        #topCoreg = createTOPSSwathSLCProduct()
        topCoreg = ut.coregSwathSLCProduct()
        topCoreg.configure()

        if inps.overlap:    
            botCoreg = ut.coregSwathSLCProduct()
            botCoreg.configure()

        for ii in range(minBurst, maxBurst):
            jj = secondaryBurstStart + ii - minBurst

            if inps.overlap:
                botBurst = referenceBottom.bursts[ii]
                topBurst = referenceTop.bursts[ii]
            else:
                topBurst = reference.bursts[ii]


            slvBurst = secondary.bursts[jj]

            #####Top burst processing
            try:
                offset = relShifts[jj]
            except:
                raise Exception('Trying to access shift for secondary burst index {0}, which may not overlap with reference'.format(jj))
        
            if inps.overlap:         
                outname = os.path.join(outdir, 'burst_top_%02d_%02d.slc'%(ii+1,ii+2))
        
                ####Setup initial polynomials
                ### If no misregs are given, these are zero
                ### If provided, can be used for resampling without running to geo2rdr again for fast results
                rdict = {'azpoly' : apoly,
                     'rgpoly' : rpoly,
                     'rangeOff' : os.path.join(offdir, 'range_top_%02d_%02d.off'%(ii+1,ii+2)),
                     'azimuthOff': os.path.join(offdir, 'azimuth_top_%02d_%02d.off'%(ii+1,ii+2))}

        
            ###For future - should account for azimuth and range misreg here .. ignoring for now. 
                azCarrPoly, dpoly = secondary.estimateAzimuthCarrierPolynomials(slvBurst, offset = -1.0 * offset)
        
                rdict['carrPoly'] = azCarrPoly
                rdict['doppPoly'] = dpoly
        
                outimg = resampSecondary(topBurst, slvBurst, rdict, outname, (not inps.noflat))
        
                copyBurst = copy.deepcopy(topBurst)
                ut.adjustValidSampleLine(copyBurst)
                copyBurst.image.filename = outimg.filename 
                print('After: ', copyBurst.firstValidLine, copyBurst.numValidLines)
                topCoreg.bursts.append(copyBurst)
            #######################################################

        
                slvBurst = secondary.bursts[jj+1]
                outname = os.path.join(outdir, 'burst_bot_%02d_%02d.slc'%(ii+1,ii+2))
        
        ####Setup initial polynomials
        ### If no misregs are given, these are zero
        ### If provided, can be used for resampling without running to geo2rdr again for fast results
                rdict = {'azpoly' : apoly,
                     'rgpoly' : rpoly,
                     'rangeOff' : os.path.join(offdir, 'range_bot_%02d_%02d.off'%(ii+1,ii+2)),
                     'azimuthOff': os.path.join(offdir, 'azimuth_bot_%02d_%02d.off'%(ii+1,ii+2))}
        
                azCarrPoly, dpoly = secondary.estimateAzimuthCarrierPolynomials(slvBurst, offset = -1.0 * offset)
        
                rdict['carrPoly'] = azCarrPoly
                rdict['doppPoly'] = dpoly
        
                outimg = resampSecondary(botBurst, slvBurst, rdict, outname, (not inps.noflat))
        
                copyBurst = copy.deepcopy(botBurst)
                ut.adjustValidSampleLine(copyBurst)
                copyBurst.image.filename = outimg.filename
                print('After: ', copyBurst.firstValidLine, copyBurst.numValidLines)
                botCoreg.bursts.append(copyBurst)
           #######################################################

            else:

                outname = os.path.join(outdir, 'burst_%02d.slc'%(ii+1))  
        
        ####Setup initial polynomials
        ### If no misregs are given, these are zero
        ### If provided, can be used for resampling without running to geo2rdr again for fast results
                rdict = {'azpoly' : apoly,
                     'rgpoly' : rpoly,
                      'rangeOff' : os.path.join(offdir, 'range_%02d.off'%(ii+1)),
                      'azimuthOff': os.path.join(offdir, 'azimuth_%02d.off'%(ii+1))}
                 

        ###For future - should account for azimuth and range misreg here .. ignoring for now.
                azCarrPoly, dpoly = secondary.estimateAzimuthCarrierPolynomials(slvBurst, offset = -1.0 * offset)
        
                rdict['carrPoly'] = azCarrPoly
                rdict['doppPoly'] = dpoly
        
                outimg = resampSecondary(topBurst, slvBurst, rdict, outname, (not inps.noflat))
                minAz, maxAz, minRg, maxRg = ut.getValidLines(slvBurst, rdict, outname,
                    misreg_az = misreg_az - offset, misreg_rng = misreg_rg)
                

                copyBurst = copy.deepcopy(topBurst)
                ut.adjustValidSampleLine_V2(copyBurst, slvBurst, minAz=minAz, maxAz=maxAz,
                                         minRng=minRg, maxRng=maxRg)
                copyBurst.image.filename = outimg.filename
                print('After: ', copyBurst.firstValidLine, copyBurst.numValidLines)
                topCoreg.bursts.append(copyBurst)
        #######################################################       

 
        topCoreg.numberOfBursts = len(topCoreg.bursts)
        topCoreg.source = ut.asBaseClass(secondary)

        if inps.overlap:
            botCoreg.numberOfBursts = len(botCoreg.bursts)
            topCoreg.reference = ut.asBaseClass(referenceTop)
            botCoreg.reference = ut.asBaseClass(referenceBottom)
            botCoreg.source = ut.asBaseClass(secondary)
            ut.saveProduct(topCoreg, outdir + '_top.xml')
            ut.saveProduct(botCoreg, outdir + '_bottom.xml')

        else:
            topCoreg.reference = reference
            ut.saveProduct(topCoreg, outdir + '.xml')    

if __name__ == '__main__':
    '''
    Main driver.
    '''
    # Main Driver
    main()



