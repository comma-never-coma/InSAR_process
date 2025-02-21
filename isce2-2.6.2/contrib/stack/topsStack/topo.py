#!/usr/bin/env python3

import numpy as np
import argparse
import os
import isce
import isceobj
import datetime
import sys
import s1a_isce_utils as ut
from isceobj.Planet.Planet import Planet
from zerodop.topozero import createTopozero
import multiprocessing as mp


def createParser():
    parser = argparse.ArgumentParser( description='Generates lat/lon/h and los for each pixel')
    parser.add_argument('-m', '--reference', type=str, dest='reference', required=True,
            help='Directory with the reference image')
    parser.add_argument('-d', '--dem', type=str, dest='dem', required=True,
            help='DEM to use for coregistration')
    parser.add_argument('-g', '--geom_referenceDir', type=str, dest='geom_referenceDir', default='geom_reference',
            help='Directory for geometry files of the reference')
    parser.add_argument('-n', '--numProcess', type=int, dest='numProcess', default=1,
            help='Number of parallel processes (default: %(default)s).')
    parser.add_argument('-useGPU', '--useGPU', dest='useGPU',action='store_true', default=False,
            help='Allow App to use GPU when available')

    return parser

def cmdLineParse(iargs=None):
    '''
    Command line parser.
    '''
    parser = createParser()
    return parser.parse_args(args=iargs)


def call_topo(input):

    (dirname, demImage, reference, ind) = input

    burst = reference.bursts[ind]
    latname = os.path.join(dirname, 'lat_%02d.rdr' % (ind + 1))
    lonname = os.path.join(dirname, 'lon_%02d.rdr' % (ind + 1))
    hgtname = os.path.join(dirname, 'hgt_%02d.rdr' % (ind + 1))
    losname = os.path.join(dirname, 'los_%02d.rdr' % (ind + 1))
    maskname = os.path.join(dirname, 'shadowMask_%02d.rdr' % (ind + 1))
    incname = os.path.join(dirname, 'incLocal_%02d.rdr' % (ind + 1))
    #####Run Topo
    planet = Planet(pname='Earth')
    topo = createTopozero()
    topo.slantRangePixelSpacing = burst.rangePixelSize
    topo.prf = 1.0 / burst.azimuthTimeInterval
    topo.radarWavelength = burst.radarWavelength
    topo.orbit = burst.orbit
    topo.width = burst.numberOfSamples
    topo.length = burst.numberOfLines
    topo.wireInputPort(name='dem', object=demImage)
    topo.wireInputPort(name='planet', object=planet)
    topo.numberRangeLooks = 1
    topo.numberAzimuthLooks = 1
    topo.lookSide = -1
    topo.sensingStart = burst.sensingStart
    topo.rangeFirstSample = burst.startingRange
    topo.demInterpolationMethod = 'BIQUINTIC'
    topo.latFilename = latname
    topo.lonFilename = lonname
    topo.heightFilename = hgtname
    topo.losFilename = losname
    topo.maskFilename = maskname
    topo.incFilename = incname
    topo.topo()

    bbox = [topo.minimumLatitude, topo.maximumLatitude, topo.minimumLongitude, topo.maximumLongitude]

    topo = None

    return bbox


def main(iargs=None):

    inps = cmdLineParse(iargs)
    
    swathList = ut.getSwathList(inps.reference)

    # see if the user compiled isce with GPU enabled
    run_GPU = False
    try:
        from zerodop.GPUtopozero.GPUtopozero import PyTopozero
        run_GPU = True
    except:
        pass
        
    if inps.useGPU and not run_GPU:
        print("GPU mode requested but no GPU ISCE code found")
    
    # setting the respective version of geo2rdr for CPU and GPU
    if run_GPU and inps.useGPU:
        print('GPU mode')
        
        from zerodop.GPUtopozero.GPUtopozero import PyTopozero
        from isceobj import Constants as CN
        from isceobj.Util.Poly2D import Poly2D
        from iscesys import DateTimeUtil as DTU
        
        frames = []
        swaths = []
        
        for swath in swathList:
            reference =  ut.loadProduct(os.path.join(inps.reference , 'IW{0}.xml'.format(swath)))
            frames.append(reference)
            swaths.append(swath)
                
        
        if len(frames) == 0:
            raise Exception('There is no common region between the two dates to process')
        
        topSwath = min(frames, key=lambda x: x.sensingStart)
        leftSwath = min(frames, key=lambda x: x.startingRange)
        bottomSwath = max(frames, key=lambda x: x.sensingStop)
        rightSwath = max(frames, key=lambda x: x.farRange)
        
        r0 = leftSwath.startingRange
        rmax = rightSwath.farRange
        dr = frames[0].bursts[0].rangePixelSize
        t0 = topSwath.sensingStart
        tmax = bottomSwath.sensingStop
        dt = frames[0].bursts[0].azimuthTimeInterval
        wvl = frames[0].bursts[0].radarWavelength
        width = int(np.round((rmax-r0)/dr) + 1)
        lgth = int(np.round((tmax-t0).total_seconds()/dt) + 1)
        
        polyDoppler = Poly2D(name='topsStack_dopplerPoly')
        polyDoppler.setWidth(width)
        polyDoppler.setLength(lgth)
        polyDoppler.setNormRange(1.0)
        polyDoppler.setNormAzimuth(1.0)
        polyDoppler.setMeanRange(0.0)
        polyDoppler.setMeanAzimuth(0.0)
        polyDoppler.initPoly(rangeOrder=0,azimuthOrder=0, coeffs=[[0.]])
        polyDoppler.createPoly2D()
        
        slantRangeImage = Poly2D()
        slantRangeImage.setWidth(width)
        slantRangeImage.setLength(lgth)
        slantRangeImage.setNormRange(1.0)
        slantRangeImage.setNormAzimuth(1.0)
        slantRangeImage.setMeanRange(0.)
        slantRangeImage.setMeanAzimuth(0.)
        slantRangeImage.initPoly(rangeOrder=1,azimuthOrder=0,coeffs=[[r0,dr]])
        slantRangeImage.createPoly2D()
        
        dirname = inps.geom_referenceDir
        os.makedirs(dirname, exist_ok=True)
        
        latImage = isceobj.createImage()
        latImage.initImage(os.path.join(dirname, 'lat.rdr'), 'write', width, 'DOUBLE')
        latImage.createImage()

        lonImage = isceobj.createImage()
        lonImage.initImage(os.path.join(dirname, 'lon.rdr'), 'write', width, 'DOUBLE')
        lonImage.createImage()

        losImage = isceobj.createImage()
        losImage.initImage(os.path.join(dirname, 'los.rdr'), 'write', width, 'FLOAT', bands=2, scheme='BIL')
        losImage.setCaster('write', 'DOUBLE')
        losImage.createImage()

        heightImage = isceobj.createImage()
        heightImage.initImage(os.path.join(dirname, 'hgt.rdr'),'write', width, 'DOUBLE')
        heightImage.createImage()
        
        incImage = isceobj.createImage()
        incImage.initImage(os.path.join(dirname, 'inc.rdr'), 'write', width, 'FLOAT', bands=2, scheme='BIL')
        incImage.setCaster('write', 'DOUBLE')
        incImage.createImage()

        #maskImage = isceobj.createImage()
        #maskImage.initImage(os.path.join(dirname, 'mask.rdr'), 'write', width, 'BYTE', bands=1, scheme='BIL')
        #maskImage.setCaster('write', 'DOUBLE')
        #maskImage.createImage()

        demImage = isceobj.createDemImage()
        demImage.load(inps.dem + '.xml')
        demImage.setCaster('read', 'FLOAT')
        demImage.createImage()
        
        orb = getMergedOrbit(frames)
        pegHdg = np.radians(orb.getENUHeading(t0))
        
        elp = Planet(pname='Earth').ellipsoid
        
        topo = PyTopozero()
        topo.set_firstlat(demImage.getFirstLatitude())
        topo.set_firstlon(demImage.getFirstLongitude())
        topo.set_deltalat(demImage.getDeltaLatitude())
        topo.set_deltalon(demImage.getDeltaLongitude())
        topo.set_major(elp.a)
        topo.set_eccentricitySquared(elp.e2)
        topo.set_rSpace(dr)
        topo.set_r0(r0)
        topo.set_pegHdg(pegHdg)
        topo.set_prf(1.0/dt)
        topo.set_t0(DTU.seconds_since_midnight(t0))
        topo.set_wvl(wvl)
        topo.set_thresh(.05)
        topo.set_demAccessor(demImage.getImagePointer())
        topo.set_dopAccessor(polyDoppler.getPointer())
        topo.set_slrngAccessor(slantRangeImage.getPointer())
        topo.set_latAccessor(latImage.getImagePointer())
        topo.set_lonAccessor(lonImage.getImagePointer())
        topo.set_losAccessor(losImage.getImagePointer())
        topo.set_heightAccessor(heightImage.getImagePointer())
        topo.set_incAccessor(incImage.getImagePointer())
        topo.set_maskAccessor(0)
        topo.set_numIter(25)
        topo.set_idemWidth(demImage.getWidth())
        topo.set_idemLength(demImage.getLength())
        topo.set_ilrl(-1)
        topo.set_extraIter(10)
        topo.set_length(lgth)
        topo.set_width(width)
        topo.set_nRngLooks(1)
        topo.set_nAzLooks(1)
        topo.set_demMethod(5) # BIQUINTIC METHOD
        topo.set_orbitMethod(0) # HERMITE
        
        # Need to simplify orbit stuff later
        nvecs = len(orb._stateVectors)
        topo.set_orbitNvecs(nvecs)
        topo.set_orbitBasis(1) # Is this ever different?
        topo.createOrbit() # Initializes the empty orbit to the right allocated size
        count = 0
        for sv in orb._stateVectors:
            td = DTU.seconds_since_midnight(sv.getTime())
            pos = sv.getPosition()
            vel = sv.getVelocity()
            topo.set_orbitVector(count,td,pos[0],pos[1],pos[2],vel[0],vel[1],vel[2])
            count += 1

        topo.runTopo()
        
        latImage.addDescription('Pixel-by-pixel latitude in degrees.')
        latImage.finalizeImage()
        latImage.renderHdr()

        lonImage.addDescription('Pixel-by-pixel longitude in degrees.')
        lonImage.finalizeImage()
        lonImage.renderHdr()

        heightImage.addDescription('Pixel-by-pixel height in meters.')
        heightImage.finalizeImage()
        heightImage.renderHdr()
        
        descr = '''Two channel Line-Of-Sight geometry image (all angles in degrees). Represents vector drawn from target to platform.
                Channel 1: Incidence angle measured from vertical at target (always +ve).
                Channel 2: Azimuth angle measured from North in Anti-clockwise direction.'''
        losImage.setImageType('bil')
        losImage.addDescription(descr)
        losImage.finalizeImage()
        losImage.renderHdr()

        demImage.finalizeImage()
        
        descr = '''Two channel angle file.
                Channel 1: Angle between ray to target and the vertical at the sensor
                Channel 2: Local incidence angle accounting for DEM slope at target'''
        incImage.addDescription(descr)
        incImage.finalizeImage()
        incImage.renderHdr()

        #descr = 'Radar shadow-layover mask. 1 - Radar Shadow. 2 - Radar Layover. 3 - Both.'
        #maskImage.addDescription(descr)
        #maskImage.finalizeImage()
        #maskImage.renderHdr()
        
        
        if slantRangeImage:
            try:
                slantRangeImage.finalizeImage()
            except:
                pass
                
                
        ####Start creating VRTs to point to global topo output
        for swath, frame in zip(swaths, frames):
            outname = os.path.join(dirname, 'IW{0}'.format(swath))
            os.makedirs(outname, exist_ok=True)
            
            for ind, burst in enumerate(frame.bursts):
                top = int(np.rint((burst.sensingStart - t0).total_seconds()/dt))
                bottom = top + burst.numberOfLines
                left = int(np.rint((burst.startingRange - r0)/dr))
                right = left + burst.numberOfSamples
                
                buildVRT( os.path.join(dirname, 'lat.rdr'),
                          os.path.join(outname, 'lat_%02d.rdr'%(ind+1)),
                          [width, lgth],
                          [top,bottom, left, right],
                          bands=1,
                          dtype='DOUBLE')

                buildVRT( os.path.join(dirname, 'lon.rdr'),
                          os.path.join(outname, 'lon_%02d.rdr'%(ind+1)),
                          [width, lgth],
                          [top,bottom, left, right],
                          bands=1,
                          dtype='DOUBLE')

                buildVRT( os.path.join(dirname, 'hgt.rdr'),
                          os.path.join(outname, 'hgt_%02d.rdr'%(ind+1)),
                          [width, lgth],
                          [top,bottom, left, right],
                          bands=1,
                          dtype='DOUBLE')

                buildVRT( os.path.join(dirname, 'los.rdr'),
                          os.path.join(outname, 'los_%02d.rdr'%(ind+1)),
                          [width, lgth],
                          [top,bottom, left, right],
                          bands=2,
                          dtype='FLOAT')
                          
                buildVRT( os.path.join(dirname, 'inc.rdr'),
                          os.path.join(outname, 'incLocal_%02d.rdr'%(ind+1)),
                          [width, lgth],
                          [top,bottom, left, right],
                          bands=2,
                          dtype='FLOAT')
            

    else:
        print('CPU mode')
        demImage = isceobj.createDemImage()
        demImage.load(inps.dem + '.xml')
        boxes = []
        inputs = []

        for swath in swathList:
            reference =  ut.loadProduct(os.path.join(inps.reference , 'IW{0}.xml'.format(swath)))
            
            ###Check if geometry directory already exists.
            dirname = os.path.join(inps.geom_referenceDir, 'IW{0}'.format(swath))
            os.makedirs(dirname, exist_ok=True)

            for ind in range(reference.numberOfBursts):
                inputs.append((dirname, demImage, reference, ind))

        # parallel processing
        print('running in parallel with {} processes'.format(inps.numProcess))
        pool = mp.Pool(inps.numProcess)
        results = pool.map(call_topo, inputs)
        pool.close()

        for bbox in results:
            boxes.append(bbox)

        boxes = np.array(boxes)
        bbox = [np.min(boxes[:,0]), np.max(boxes[:,1]), np.min(boxes[:,2]), np.max(boxes[:,3])]
        print('bbox : ', bbox)
    

if __name__ == '__main__':

    main()

def getMergedOrbit(product):
        from isceobj.Orbit.Orbit import Orbit

        ###Create merged orbit
        orb = Orbit()
        orb.configure()

        burst = product[0].bursts[0]
        #Add first burst orbit to begin with
        for sv in burst.orbit:
             orb.addStateVector(sv)


        for pp in product:
            ##Add all state vectors
            for bb in pp.bursts:
                for sv in bb.orbit:
                    if (sv.time< orb.minTime) or (sv.time > orb.maxTime):
                        orb.addStateVector(sv)

                bb.orbit = orb

        return orb
        
def buildVRT(srcname, dstname, dims, bbox, bands=1, dtype='FLOAT'):
    '''
    Write a VRT to point to the parent mosaicked file.
    '''

    header='<VRTDataset rasterXSize="{width}" rasterYSize="{lgth}">'
    band = '''    <VRTRasterBand dataType="{dtype}" band="{band}">
        <NoDataValue>0.0</NoDataValue>
        <SimpleSource>
            <SourceFilename relativeToVRT="1">{relpath}</SourceFilename>
            <SourceBand>{band}</SourceBand>
            <SourceProperties RasterXSize="{gwidth}" RasterYSize="{glgth}" DataType="{dtype}"/>
            <SrcRect xOff="{left}" yOff="{top}" xSize="{width}" ySize="{lgth}"/>
            <DstRect xOff="0" yOff="0" xSize="{width}" ySize="{lgth}"/>
        </SimpleSource>
    </VRTRasterBand>
'''
    tail = "</VRTDataset>"

    width = bbox[3] - bbox[2]
    lgth = bbox[1] - bbox[0]

    odtype = dtype
    if dtype.upper() == 'FLOAT':
        dtype = 'Float32'
    elif dtype.upper() == 'DOUBLE':
        dtype = 'Float64'
    elif dtype.upper() == 'BYTE':
        dtype = 'UInt8'
    else:
        raise Exception('Unsupported type {0}'.format(dtype))

    relpath = os.path.relpath(srcname + '.vrt', os.path.dirname(dstname))
    gwidth = dims[0]
    glgth = dims[1]
    left = bbox[2]
    top = bbox[0]
    
    img = isceobj.createImage()
    img.bands = bands
    img.scheme = 'BIL'
    img.setWidth(width)
    img.setLength(lgth)
    img.dataType = odtype
    img.filename = dstname
    img.setAccessMode('READ')
    img.renderHdr()

    with open(dstname + '.vrt', 'w') as fid:
        fid.write( header.format(width=width, lgth=lgth) + '\n')

        for bnd in range(bands):
            fid.write( band.format(width=width, lgth=lgth,
                                   gwidth=gwidth, glgth=glgth,
                                   left=left, top=top,
                                   relpath=relpath, dtype=dtype,
                                   band=bnd+1))

        fid.write(tail + '\n')
