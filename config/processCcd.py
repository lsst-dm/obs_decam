from lsst.obs.mosaic.isr import MosaicIsrTask
config.isr.retarget(MosaicIsrTask)
config.charImage.repair.cosmicray.nCrPixelMax = 100000
