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
def setWcs(dataId, repoPath, exp):
    path = repoPath + "/preprocessed/%s/%s/%s/%s/%s.fits"%(dataId['field'], dataId['subfield'],
         dataId['filter'], dataId['dateObs'], dataId['objname'])
    hdus = pyfits.open(path)
    #  the ccdnum is 1 based, but so are the hdus because hdus[0] has no data
    header = hdus[dataId['ccdnum']].header
    hdus.close()
    crval1 = float(header['CRVAL1'])
    crval2 = float(header['CRVAL2'])
    crpix1 = float(header['CRPIX1'])
    crpix2 = float(header['CRPIX2'])
    cd11 = float(header['CD1_1'])
    cd21 = float(header['CD2_1'])
    cd12 = float(header['CD1_2'])
    cd22 = float(header['CD2_2'])
    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']
    crval = lsst.afw.coord.Coord(lsst.afw.geom.PointD(crval1,crval2)) #((139.04807,30.03947)
    crpix = lsst.afw.geom.PointD(crpix1,crpix2) #(-44.0659, 4107.8164)
    wcs = lsst.afw.image.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)
    xc = (naxis1/2)-1
    yc = (naxis2/2)-1
    radec = wcs.pixelToSky(xc,yc)
    ra = radec.toFk5().getRa().asDegrees()
    dec = radec.toFk5().getDec().asDegrees()
    crval = lsst.afw.coord.Coord(lsst.afw.geom.PointD(ra,dec)) #((139.04807,30.03947)
    crpix = lsst.afw.geom.PointD(xc,yc) #(-44.0659, 4107.8164)
    wcs = lsst.afw.image.makeWcs(crval,crpix,cd11,cd12,cd21,cd22)
    exp.setWcs(wcs)

#   Set the mask plane for the "BAD" pixels using the BPM supplied by the observatory
def setMask(dataId, repoPath, exp):
    path = repoPath + "/preprocessed/%s/%s/%s/%s/%s.fits"%(dataId['field'], dataId['subfield'],
         dataId['filter'], dataId['dateObs'], dataId['objname'])
    hdus = pyfits.open(path)
    #  the ccdnum is 1 based, but so are the hdus because hdus[0] has no data
    header = hdus[0].header
    observat = header['OBSERVAT'].lower()
    hdus.close()
    mask = exp.getMaskedImage().getMask()
    #   First set the saturated pixels with the SAT bit
    satbitm = mask.getPlaneBitMask('SAT')
    satmask = (exp.getMaskedImage().getImage().getArray() >= 20000)
    maskarray = mask.getArray()
    mask.getArray()[satmask] |= satbitm
    #   First get the bpm for that date, ccd, and telescope
    bpmpath = "/home/DLS/masks/%s/%s/bpm_%d.fits"%(observat, dataId["dateObs"][0:4],
                dataId['ccdnum'])
    if os.path.isfile(bpmpath):
        badbitm = mask.getPlaneBitMask('BAD')
        img = lsst.afw.image.ImageI(bpmpath)
        badmask = (img.getArray() < 1)
        mask.getArray()[badmask] |= badbitm
    #    Next get the exclusion regions for each exposure, mark as "SUSPECT"
    regpath = "/home/DLS/masks/%s/%s/%s/%s/%s_%d.bpm.fits.fz"%(dataId['field'], dataId['subfield'],
         dataId['filter'], dataId['dateObs'], dataId['objname'], dataId['ccdnum'])
    if os.path.isfile(regpath):
        badbitm = mask.getPlaneBitMask('SUSPECT')
        regimg = lsst.afw.image.ImageI(regpath)
        badmask = (regimg.getArray() < 1)
        mask.getArray()[badmask] |= badbitm


def setupButler(root):
    butler=Butler(root=root)
    return butler

def updateVar(exp):
    ccd=exp.getDetector()[0]
    gain=ccd.getGain()
    readNoise=ccd.getReadNoise()
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
    exp = butler.get('preprocessed', dataId)
    setWcs(dataId, repoPath, exp)
    setMask(dataId, repoPath, exp)
    updateVar(exp)
    butler.put(exp, 'postISRCCD', dataId)
    charIm(butler, exp, dataId)

if __name__ == '__main__':
    run()
