import os.path
from lsst.utils import getPackageDir
from lsst.meas.extensions.astrometryNet import LoadAstrometryNetObjectsTask

config.match.refObjLoader.retarget(LoadAstrometryNetObjectsTask)
config.match.refObjLoader.load(os.path.join(getPackageDir("obs_mosaic"), "config", "filterMap.py"))
config.match.refObjLoader.defaultFilter = 'r'
