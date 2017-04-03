from builtins import map
#
# LSST Data Management System
# Copyright 2012-2016 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import os
import re
import numpy as np
import pyfits
from lsst.utils import getPackageDir
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
from lsst.obs.base import CameraMapper, exposureFromImage
from lsst.daf.persistence import ButlerLocation
from lsst.ip.isr import isr
import lsst.pex.policy as pexPolicy
from .makeMosaicRawVisitInfo import MakeMosaicRawVisitInfo

np.seterr(divide="ignore")

class MosaicMapper(CameraMapper):
    packageName = 'obs_mosaic'

    MakeRawVisitInfoClass = MakeMosaicRawVisitInfo

    detectorNames = {1:'E1', 2:'E2', 3:'E3', 4:'E4', 5:'W5', 6:'W6', 7:'W7', 8:'W8'}

    def __init__(self, inputPolicy=None, **kwargs):
        policyFile = pexPolicy.DefaultPolicyFile(self.packageName, "MosaicMapper.paf", "policy")
        policy = pexPolicy.Policy(policyFile)

        super(MosaicMapper, self).__init__(policy, policyFile.getRepositoryPath(), **kwargs)

        #I found these values in the mosaic 1 manual from september 2004. lambda is in nm
        afwImageUtils.defineFilter('B', lambdaEff=436, alias=['B'])
        afwImageUtils.defineFilter('V', lambdaEff=537, alias=['V'])
        afwImageUtils.defineFilter('R', lambdaEff=644, alias=['R'])
        afwImageUtils.defineFilter('z', lambdaEff=940, alias=["SDSS z'"])

        # The data ID key ccdnum is not directly used in the current policy
        # template of the raw dataset, so is not in its keyDict automatically.
        # Add it so raw dataset know about the data ID key ccdnum.
        self.mappings["raw"].keyDict.update({'ccdnum': int})

        # The number of bits allocated for fields in object IDs
        # TODO: This needs to be updated; also see Trac #2797
        MosaicMapper._nbit_tract = 10
        MosaicMapper._nbit_patch = 10
        MosaicMapper._nbit_filter = 4
        MosaicMapper._nbit_id = 64 - (MosaicMapper._nbit_tract +
                                     2*MosaicMapper._nbit_patch +
                                     MosaicMapper._nbit_filter)

    def _extractDetectorName(self, dataId):
        copyId = self._transformId(dataId)
        try:
            return MosaicMapper.detectorNames[copyId['ccdnum']]
        except KeyError:
            raise RuntimeError("No name found for dataId: %s"%(dataId))

    def _transformId(self, dataId):
        copyId = CameraMapper._transformId(self, dataId)
        if "ccd" in copyId:
            copyId.setdefault("ccdnum", copyId["ccd"])
        return copyId

    def bypass_ccdExposureId(self, datasetType, pythonType, location, dataId):
        return self._computeCcdExposureId(dataId)

    def bypass_ccdExposureId_bits(self, datasetType, pythonType, location, dataId):
        return 32  # not really, but this leaves plenty of space for sources

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.

        @param dataId (dict) Data identifier with visit, ccd
        """
        copyId = self._transformId(dataId)
        date=dataId['dateObs'].replace('-','')
        #you have to do remove the middle 2 didgets in the date so the
        #exposureId ddoesnt exceed 32 bits and is still unique for each ccd exp
        dateCat=date[0]+date[3]+date[4:]
        obj=dataId['objname'].replace('obj','')
        ccdnum = copyId['ccdnum']
        return int("%03s%01s" % (obj, ccdnum))

    def _computeCoaddExposureId(self, dataId, singleFilter):
        """Compute the 64-bit (long) identifier for a coadd.

        @param dataId (dict)       Data identifier with tract and patch.
        @param singleFilter (bool) True means the desired ID is for a single-
                                   filter coadd, in which case dataId
                                   must contain filter.
        """
        tract = int(dataId['tract'])
        if tract < 0 or tract >= 2**MosaicMapper._nbit_tract:
            raise RuntimeError('tract not in range [0,%d)' % (2**MosaicMapper._nbit_tract))
        patchX, patchY = [int(x) for x in dataId['patch'].split(',')]
        for p in (patchX, patchY):
            if p < 0 or p >= 2**MosaicMapper._nbit_patch:
                raise RuntimeError('patch component not in range [0, %d)' % 2**MosaicMapper._nbit_patch)
        oid = (((tract << MosaicMapper._nbit_patch) + patchX) << MosaicMapper._nbit_patch) + patchY
        if singleFilter:
            return (oid << MosaicMapper._nbit_filter) + afwImage.Filter(dataId['filter']).getId()
        return oid

    def bypass_deepCoaddId(self, datasetType, pythonType, location, dataId):
        return self._computeCoaddExposureId(dataId, True)

    def bypass_deepCoaddId_bits(self, *args, **kwargs):
        return 64 - MosaicMapper._nbit_id

    def bypass_deepMergedCoaddId(self, datasetType, pythonType, location, dataId):
        return self._computeCoaddExposureId(dataId, False)

    def bypass_deepMergedCoaddId_bits(self, *args, **kwargs):
        return 64 - MosaicMapper._nbit_id

    def std_preprocessed(self, item, dataId):
        """Standardize a preprocess dataset by converting it to an Exposure.

        Preprocessed images are MEF files with one HDU for each detector.
        Header keyword MJD-OBS exist only in the zeroth
        extension and is copied to metadata.

        @param item: The image read by the butler
        @param dataId: Data identifier
        @return (lsst.afw.image.Exposure) the standardized Exposure
        """
        # Convert the raw DecoratedImage to an Exposure, set metadata and wcs.
        md = item.getMetadata()

        # Convert wcs to a TAN by removing the higher order corrections which
        # were put by the DLS pipeline into the image Wcs.  The Wcs which results
        # is now only approximate
        md.set('CTYPE1', 'RA---TAN')
        md.set('CTYPE2', 'DEC--TAN')
        for kw in md.paramNames():
            if kw.startswith('WAT1') or kw.startswith('WAT2'):
                md.remove(kw)
        exp = exposureFromImage(item)

        #   convert the hdu0 header to visitInfo using pyfits to read it
        path = self.map_preprocessed(dataId).getLocationsWithRoot()[0]
        headerPath = re.sub(r'[\[](\d+)[\]]$', "", path)
        header = pyfits.open(headerPath)[0].header
        #  We don't actually use all of thes values in the header, but put this
        #  list here for later reference
        md0 = type(md)()
        extractKeys = ('DATE', 'FILENAME', 'EXPTIME', 'DARKTIME', 'RA', 'DEC', 'DATE-OBS',
                       'TIME-OBS', 'MJD-OBS', 'OBSERVAT', 'TELESCOP', 'TELRADEC', 'TELRA',
                       'TELDEC', 'ZD', 'AIRMASS', 'DETECTOR', 'FILTER', 'READTIME', 'OBSID')
        for key in extractKeys:
            if key in header.keys():
                md0.add(key, header[key])
        #   TIMESYS is utc approximate in the header, so we need to replace it
        md0.add('TIMESYS', 'utc')
        exposureId = self._computeCcdExposureId(dataId)
        visitInfo = self.makeRawVisitInfo(md=md0, exposureId=exposureId)
        exp.getInfo().setVisitInfo(visitInfo)
        # Keywords EXPTIME and MJD-OBS are used to set the calib object.
        for kw in ('MJD-OBS', 'EXPTIME', 'OBSERVAT'):
            if kw in md0.paramNames():
                md.add(kw, md0.get(kw))
        # Standardize an Exposure, including setting the calib object
        result = self._standardizeExposure(self.exposures['preprocessed'], exp, dataId,
                                         trimmed=False)
        return result
