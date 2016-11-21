import pyfits
import sys
from lsst.utils import getPackageDir
from lsst.daf.persistence import Butler
from lsst.ip.isr.isr import updateVariance, makeThresholdMask, maskPixelsFromDefectList
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask
import lsst.afw.coord
import lsst.afw.geom
import lsst.afw.image
import os
import pdb

#  Use the header from the preprocessed mosaic image to set the wcs of the exposure.
#  The wcs is centered on the central pixel, using the coordinate
#  transform in the header of the image for the specified CCD.
def makeWcs(metadata, dimensions):
    crval1 = float(metadata.get('CRVAL1'))
    crval2 = float(metadata.get('CRVAL2'))
    crpix1 = float(metadata.get('CRPIX1'))
    crpix2 = float(metadata.get('CRPIX2'))
    cd11 = float(metadata.get('CD1_1'))
    cd21 = float(metadata.get('CD2_1'))
    cd12 = float(metadata.get('CD1_2'))
    cd22 = float(metadata.get('CD2_2'))
    crval = lsst.afw.coord.Coord(lsst.afw.geom.PointD(crval1,crval2)) #((139.04807,30.03947)
    crpix = lsst.afw.geom.PointD(crpix1,crpix2) #(-44.0659, 4107.8164)
    wcs = lsst.afw.image.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)
    xc = dimensions.getX()
    yc = dimensions.getY()
    radec = wcs.pixelToSky(xc,yc)
    ra = radec.toFk5().getRa().asDegrees()
    dec = radec.toFk5().getDec().asDegrees()
    crval = lsst.afw.coord.Coord(lsst.afw.geom.PointD(ra,dec)) #((139.04807,30.03947)
    crpix = lsst.afw.geom.PointD(xc,yc) #(-44.0659, 4107.8164)
    wcs = lsst.afw.image.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)
    return wcs

#   Set the mask plane for the "BAD" pixels using the BPM supplied by the observatory
def setMask(butler, dataId, exp):
    mask = exp.getMaskedImage().getMask()
    #   First set the saturated pixels with the SAT bit
    satbitm = mask.getPlaneBitMask('SAT')
    satmask = (exp.getMaskedImage().getImage().getArray() >= 20000)
    maskarray = mask.getArray()
    mask.getArray()[satmask] |= satbitm

    #   First get the bpm for that date, ccd, and telescope
    dataId['observatory'] = 'kpno'
    dataId['year'] = dataId['dateObs'][0:4]
    if butler.datasetExists('bpm', dataId):
        img = butler.get('bpm', dataId)
        badbitm = mask.getPlaneBitMask('BAD')
        badmask = (img.getArray() < 1)
        mask.getArray()[badmask] |= badbitm

    #    Next get the exclusion regions for each exposure, mark as "SUSPECT"
    regpath = "/home/DLS/masks/%s/%s/%s/%s/%s_%d.bpm.fits.fz"%(dataId['field'], dataId['subfield'],
         dataId['filter'], dataId['dateObs'], dataId['objname'], dataId['ccdnum'])
    if butler.datasetExists('masked', dataId):
        regimg = butler.get('masked', dataId)
        badbitm = mask.getPlaneBitMask('SUSPECT')
        badmask = (regimg.getArray() < 1)
        mask.getArray()[badmask] |= badbitm


def setupButler(root):
    butler=Butler(root=root)
    return butler

def updateVar(exp, metadata):
    gain=metadata.get('GAIN')
    readNoise=metadata.get('RDNOISE')
    readNoise_ADU=readNoise/gain
    updateVariance(exp.getMaskedImage(), gain, readNoise_ADU)

def charIm(butler,exp, dataId):
    dataRef=butler.dataRef('preprocessed', dataId=dataId)
    charIm=CharacterizeImageTask()
    charIm.run(dataRef, exposure=exp, doUnpersist=False)

def run():
    datadir= getPackageDir('testdata_mosaic')
    repoPath=os.path.join(datadir, 'rawData')
    butler=setupButler(repoPath)
    dataId={"field":'F2', "subfield":'p23', "filter":"R", "dateObs":'2003-01-04', "objname":'obj074', "ccdnum": 1}
    #   The preprocessed image is one plane of the mosaic preprocessed by DLS
    #   It us read as a Decorated Image so that the metadata is available
    decimg = butler.get('preprocessed', dataId)
    mask = lsst.afw.image.MaskU(decimg.getDimensions())
    var = lsst.afw.image.ImageF(decimg.getDimensions())
    exp = lsst.afw.image.ExposureF(lsst.afw.image.MaskedImageF(decimg.getImage(), mask, var))

    #   Create a WCS centered on the ccd image extent
    exp.setWcs(makeWcs(decimg.getMetadata(), decimg.getDimensions()))

    #   Use the butler to fetch info needed to use the DLS bpm and masked region files
    setMask(butler, dataId, exp)

    #   Update the variance plane using the image prior to background subtraction
    updateVar(exp, decimg.getMetadata())

    #   Save as a postISRCCD
    butler.put(exp, 'postISRCCD', dataId)
    charIm(butler, exp, dataId)

if __name__ == '__main__':
    run()
