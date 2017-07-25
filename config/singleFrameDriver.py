import os

from lsst.utils import getPackageDir

config.processCcd.load(os.path.join(getPackageDir("obs_mosaic"), "config", "processCcd.py"))
config.ccdKey = 'ccdnum'
