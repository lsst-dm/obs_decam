from __future__ import absolute_import, division, print_function
# Enable CModel mags (unsetup meas_modelfit to disable)
# 'config' is a SourceMeasurementConfig.
#this is copied from obs_subaru/config. i commented out their change to
#the classification extendedness value. i'm not sure if that is HSC specifit at the moment
import os
try:
    import lsst.meas.modelfit
    config.measurement.plugins.names |= ["modelfit_DoubleShapeletPsfApprox", "modelfit_CModel"]
    config.measurement.slots.modelFlux = 'modelfit_CModel'
#    config.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio = 0.985
except (KeyError, ImportError):
    print("Cannot import lsst.meas.modelfit: disabling CModel measurements")
