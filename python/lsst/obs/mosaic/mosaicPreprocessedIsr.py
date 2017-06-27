#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import lsst.afw.coord
import lsst.afw.geom
import lsst.afw.image
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask
from lsst.ip.isr.isrFunctions import updateVariance, makeThresholdMask, maskPixelsFromDefectList, interpolateFromMask

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
    mi = exp.getMaskedImage()
    mask = mi.getMask()
    #   First set the saturated pixels with the SAT bit
    satbitm = mask.getPlaneBitMask('SAT')
    satlevel = exp.getDetector()[0].getSaturation() * .80
    satmask = (exp.getMaskedImage().getImage().getArray() >= satlevel)
    maskarray = mask.getArray()
    mask.getArray()[satmask] |= satbitm

    #   First get the bpm for that date, ccd, and telescope
    dI2 = dict(dataId)
    dI2['observatory'] = exp.getMetadata().get('OBSERVAT').lower()
    dI2['year'] = dI2['dateObs'][0:4]
    if butler.datasetExists('bpm', dI2):
        img = butler.get('bpm', dI2)
        badbitm = mask.getPlaneBitMask('BAD')
        badmask = (img.getArray() == 0)
        mask.getArray()[badmask] |= badbitm

    #    Next get the exclusion regions for each exposure, mark as "SUSPECT"
    if butler.datasetExists('masked', dI2):
        regimg = butler.get('masked', dI2)
        badbitm = mask.getPlaneBitMask('SUSPECT')
        badmask = (regimg.getArray() == 0)
        mask.getArray()[badmask] |= badbitm

def updateVar(exp, metadata):
    gain=metadata.get('GAIN')
    readNoise=metadata.get('RDNOISE')
    readNoise_ADU=readNoise/gain
    updateVariance(exp.getMaskedImage(), gain, readNoise_ADU)

class MosaicPreprocessedIsrConfig(pexConfig.Config):
    doWrite = pexConfig.Field(
        dtype=bool,
        doc="Persist loaded data as a postISRCCD? The default is false, to avoid duplicating data.",
        default=True,
    )
    datasetType = pexConfig.Field(
        dtype=str,
        doc="Dataset type for input data; read by ProcessCcdTask; users will typically leave this alone",
        default="preprocessed",
    )

## \addtogroup LSST_task_documentation
## \{
## \page MosaicPreprocessedIsrTask
## \ref MosaicPreprocessedIsrTask_ "MosaicPreprocessedIsrTask"
## \copybrief MosaicPreprocessedIsrTask
## \}


class MosaicPreprocessedIsrTask(pipeBase.Task):
    """!Load an "instcal" exposure as a post-ISR CCD exposure

    @anchor MosaicPreprocessedIsrTask_

    @section pipe_tasks_mosaicPreprocessedIsr_Contents  Contents

     - @ref pipe_tasks_mosaicPreprocessedIsr_Purpose
     - @ref pipe_tasks_mosaicPreprocessedIsr_Initialize
     - @ref pipe_tasks_mosaicPreprocessedIsr_IO
     - @ref pipe_tasks_mosaicPreprocessedIsr_Config

    @section pipe_tasks_mosaicPreprocessedIsr_Purpose  Description

    Load "instcal" exposures from the community pipeline as a post-ISR exposure,
    and optionally persist it as a `postISRCCD`.

    This is used to retarget the `isr` subtask in `ProcessCcdTask` when you prefer to use
    the community pipeline instead of the LSST software stack to perform ISR on Mosaic images.

    @section pipe_tasks_mosaicPreprocessedIsr_Initialize  Task initialisation

    @copydoc \_\_init\_\_

    @section pipe_tasks_mosaicPreprocessedIsr_IO  Invoking the Task

    The main method is `runDataRef`.

    @section pipe_tasks_mosaicPreprocessedIsr_Config  Configuration parameters

    See @ref MosaicPreprocessedIsrConfig
    """
    ConfigClass = MosaicPreprocessedIsrConfig
    _DefaultName = "isr"

    @pipeBase.timeMethod
    def runDataRef(self, sensorRef):
        """!Load a Mosaic community pipeline "instcal" exposure as a post-ISR CCD exposure

        @param[in] sensorRef  butler data reference for post-ISR exposure
            (a daf.persistence.butlerSubset.ButlerDataRef)

        @return a pipeBase.Struct with fields:
        - exposure: exposure after application of ISR: the "instcal" exposure, unchanged
        """
        self.log.info("Loading Mosaic community pipeline file %s" % (sensorRef.dataId))
        butler = sensorRef.getButler()
        dataId = sensorRef.dataId
        exp = butler.get('preprocessed', dataId)

        #   Use the butler to fetch info needed to use the DLS bpm and masked region files
        setMask(butler, dataId, exp)
        interpolateFromMask(exp.getMaskedImage(), 1.0,  growFootprints=1, maskName='BAD')
        #   Update the variance plane using the image prior to background subtraction
        updateVar(exp, exp.getMetadata())

        #interpolateFromMask(exp.getMaskedImage(), 1.0,  growFootprints=1, maskName='SAT')
        if self.config.doWrite:
            sensorRef.put(exp, "postISRCCD")

        return pipeBase.Struct(
            exposure=exp,
        )
