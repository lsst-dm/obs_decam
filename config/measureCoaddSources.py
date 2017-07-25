"""obs mosiac specific overrides for MeasureMergedCoaddSourcesTask"""

import os.path
from lsst.utils import getPackageDir
from lsst.meas.extensions.astrometryNet import LoadAstrometryNetObjectsTask

#configs to use the astrometry refobjloader in directMatch
config.match.refObjLoader.retarget(LoadAstrometryNetObjectsTask)
config.match.refObjLoader.load(os.path.join(getPackageDir("obs_mosaic"), "config", "filterMap.py"))
config.match.refObjLoader.defaultFilter = 'r'
