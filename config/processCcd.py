from lsst.obs.mosaic.isr import MosaicIsrTask
config.isr.retarget(MosaicIsrTask)
config.charImage.repair.cosmicray.nCrPixelMax = 100000

#below overrides copied from subaru

# PSF determination
config.charImage.measurePsf.reserveFraction = 0.2
config.charImage.measurePsf.starSelector["objectSize"].sourceFluxField = 'base_PsfFlux_flux'
try:
    import lsst.meas.extensions.psfex.psfexPsfDeterminer
    config.charImage.measurePsf.psfDeterminer["psfex"].spatialOrder = 2
    config.charImage.measurePsf.psfDeterminer.name = "psfex"
except ImportError as e:
    print("WARNING: Unable to use psfex: %s" % e)
    config.charImage.measurePsf.psfDeterminer.name = "pca"

# Detection
config.charImage.detection.isotropicGrow = True
config.calibrate.detection.isotropicGrow = True

# Deblender
config.charImage.measurement.plugins.names |= ["base_Jacobian", "base_FPPosition"]
config.calibrate.measurement.plugins.names |= ["base_Jacobian", "base_FPPosition"]
