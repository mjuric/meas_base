#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
"""
Definition and registration of classification plugins
"""

import numpy

import lsst.pex.config
from .afterburner import AfterburnerPluginConfig, AfterburnerPlugin
from .pluginRegistry import register

__all__ = (
    "ABClassificationConfig", "ABClassificationPlugin",
)


class ABClassificationConfig(AfterburnerPluginConfig):

    fluxRatio = lsst.pex.config.Field(dtype=float, default=.925, optional=True,
                                      doc="critical ratio of model to psf flux")
    modelErrFactor = lsst.pex.config.Field(dtype=float, default=0.0, optional=True,
                                           doc="correction factor for modelFlux error")
    psfErrFactor = lsst.pex.config.Field(dtype=float, default=0.0, optional=True,
                                         doc="correction factor for psfFlux error")


@register("base_ClassificationExtendedness")
class ABClassificationPlugin(AfterburnerPlugin):
    """
    A binary measure of the extendedness of a source, based a simple cut on the ratio of the
    PSF flux to the model flux.

    Because the fluxes on which this algorithm is based are slot measurements, they can be provided
    by different algorithms, and the "fluxRatio" threshold used by this algorithm should generally
    be set differently for different algorithms.  To do this, plot the difference between the PSF
    magnitude and the model magnitude vs. the PSF magnitude, and look for where the cloud of galaxies
    begins.
    """

    ConfigClass = ABClassificationConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.DEFAULT_AFTERBURNER

    def __init__(self, config, name, schema, metadata):
        AfterburnerPlugin.__init__(self, config, name, schema, metadata)
        self.keyProbability = schema.addField(name + "_value", type="D",
                                              doc="Set to 1 for extended sources, 0 for point sources.")
        self.keyFlag = schema.addField(name + "_flag", type="Flag", doc="Set to 1 for any fatal failure.")

    def burn(self, measRecord):
        modelFlux = measRecord.getModelFlux()
        psfFlux = measRecord.getPsfFlux()
        modelFluxFlag = (measRecord.getModelFluxFlag()
                         if measRecord.table.getModelFluxFlagKey().isValid()
                         else False)
        psfFluxFlag = (measRecord.getPsfFluxFlag()
                       if measRecord.table.getPsfFluxFlagKey().isValid()
                       else False)
        flux1 = self.config.fluxRatio*modelFlux
        if not self.config.modelErrFactor == 0:
            flux1 += self.config.modelErrFactor*measRecord.getModelFluxErr()
        flux2 = psfFlux
        if not self.config.psfErrFactor == 0:
            flux2 += self.config.psfErrFactor*measRecord.getPsfFluxErr()

        # A generic failure occurs when either FluxFlag is set to True
        # A generic failure also occurs if either calculated flux value is NAN:
        #     this can occur if the Flux field itself is NAN,
        #     or the ErrFactor != 0 and the FluxErr is NAN
        if numpy.isnan(flux1) or numpy.isnan(flux2) or modelFluxFlag or psfFluxFlag:
            self.fail(measRecord)
        else:
            if flux1 < flux2:
                measRecord.set(self.keyProbability, 0.0)
            else:
                measRecord.set(self.keyProbability, 1.0)

    def fail(self, measRecord, error=None):
        # Override fail() to do nothing in the case of an exception.  We should be setting a flag
        # instead.
        measRecord.set(self.keyFlag, True)
