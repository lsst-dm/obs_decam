import collections
import re
from lsst.pipe.tasks.ingestCalibs import CalibsParseTask


class MosaicCalibsParseTask(CalibsParseTask):

    def getInfo(self, filename):
        """Get information about the image from the filename and/or its contents

        @param filename: Name of file to inspect
        @return File properties; list of file properties for each extension
        """
        phuInfo, infoList = CalibsParseTask.getInfo(self, filename)
        # Single-extension fits without EXTNAME can be a valid CP calibration product
        # Use info of primary header unit
        if not infoList:
            infoList.append(phuInfo)
        for info in infoList:
            info['path'] = filename
        # Try to fetch a date from filename
        # and use as the calibration dates if not already set
        found = re.search('(\d\d\d\d-\d\d-\d\d)', filename)
        if not found:
            return phuInfo, infoList
        date = found.group(1)
        for info in infoList:
            if 'calibDate' not in info or info['calibDate'] == "unknown":
                info['calibDate'] = date
        return phuInfo, infoList

    def _translateFromCalibId(self, field, md):
        """Fetch the ID from the CALIB_ID header

        Calibration products made with constructCalibs have some metadata
        saved in its FITS header CALIB_ID.
        """
        data = md.get("CALIB_ID")
        match = re.search(".*%s=(\S+)" % field, data)
        return match.groups()[0]

    def translate_ccdnum(self, md):
        """Return CCDNUM as a integer

        @param md (PropertySet) FITS header metadata
        """
        if md.exists("CCDNUM"):
            ccdnum = md.get("CCDNUM")
        else:
            return self._translateFromCalibId("ccdnum", md)
        # Some MasterCal from NOAO Archive has 2 CCDNUM keys in each HDU
        # Make sure only one integer is returned.
        if isinstance(ccdnum, collections.Sequence):
            try:
                ccdnum = ccdnum[0]
            except IndexError:
                ccdnum = None
        return ccdnum

    def translate_date(self, md):
        """Extract the date as a strong in format YYYY-MM-DD from the FITS header DATE-OBS.
        Return "unknown" if the value cannot be found or converted.

        @param md (PropertySet) FITS header metadata
        """
        if md.exists("DATE-OBS"):
            date = md.get("DATE-OBS")
            found = re.search('(\d\d\d\d-\d\d-\d\d)', date)
            if found:
                date = found.group(1)
            else:
                self.log.warn("DATE-OBS does not match format YYYY-MM-DD")
                date = "unknown"
        elif md.exists("CALIB_ID"):
            date = self._translateFromCalibId("calibDate", md)
        else:
            date = "unknown"
        return date

    def translate_filter(self, md):
        """Extract the filter name

        Translate a full filter description into a mere filter name.
        Return "unknown" if the keyword FILTER does not exist in the header,
        which can happen for some valid Community Pipeline products.

        @param md (PropertySet) FITS header metadata
        """
        if md.exists("FILTER"):
            return CalibsParseTask.translate_filter(self, md)
        elif md.exists("CALIB_ID"):
            return self._translateFromCalibId("filter", md)
        else:
            return "unknown"

    @staticmethod
    def getExtensionName(md):
        """ Get the name of the extension

        @param md (PropertySet) FITS header metadata
        """
        return md.get('EXTNAME')
