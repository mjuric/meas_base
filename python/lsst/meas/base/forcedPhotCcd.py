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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import collections

import lsst.pex.config
import lsst.pex.exceptions
import lsst.pex.logging
import lsst.pipe.base
import lsst.afw.image
import lsst.afw.table
from lsst.geom import convexHull

from .forcedPhotImage import ProcessImageForcedTask, ProcessImageForcedConfig

try:
    from lsst.meas.mosaic import applyMosaicResults
except ImportError:
    applyMosaicResults = None

__all__ = ("ForcedPhotCcdConfig", "ForcedPhotCcdTask", "PerTractCcdDataIdContainer")

class PerTractCcdDataIdContainer(lsst.pipe.base.DataIdContainer):
    """A DataIdContainer that combines raw data IDs with a tract.

    Required because we need to add "tract" to the raw data ID keys, and that's tricky.
    This IdContainer assumes that a calexp is being measured using the detection
    information from the set of coadds which intersect with the calexp. It needs the
    calexp id (e.g. visit, raft, sensor), but it also uses the tract to decide what set of
    coadds to use.  The references from the tract whose patches intersect with the calexp
    are used.
    """

    def castDataIds(self, butler):
        """Validate data IDs and cast them to the correct type (modify idList in place).

        @param butler: data butler
        """
        try:
            idKeyTypeDict = butler.getKeys(datasetType="src", level=self.level)
        except KeyError as e:
            raise KeyError("Cannot get keys for datasetType %s at level %s: %s" % ("src", self.level, e))

        idKeyTypeDict = idKeyTypeDict.copy()
        idKeyTypeDict["tract"] = int

        for dataDict in self.idList:
            for key, strVal in dataDict.iteritems():
                try:
                    keyType = idKeyTypeDict[key]
                except KeyError:
                    validKeys = sorted(idKeyTypeDict.keys())
                    raise KeyError("Unrecognized ID key %r; valid keys are: %s" % (key, validKeys))
                if keyType != str:
                    try:
                        castVal = keyType(strVal)
                    except Exception:
                        raise TypeError("Cannot cast value %r to %s for ID key %r" % (strVal, keyType, key,))
                    dataDict[key] = castVal

    def _addDataRef(self, namespace, dataId, tract):
        """Construct a dataRef based on dataId, but with an added tract key"""
        forcedDataId = dataId.copy()
        forcedDataId['tract'] = tract
        dataRef = namespace.butler.dataRef(datasetType=self.datasetType, dataId=forcedDataId)
        self.refList.append(dataRef)

    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList
        """
        if self.datasetType is None:
            raise RuntimeError("Must call setDatasetType first")
        skymap = None
        log = None
        visitTract = {} # Set of tracts for each visit
        visitRefs = {} # List of data references for each visit
        for dataId in self.idList:
            if "tract" not in dataId:
                # Discover which tracts the data overlaps
                if log is None:
                    log = lsst.pex.logging.Log.getDefaultLog()
                log.info("Reading WCS for components of dataId=%s to determine tracts" % (dict(dataId),))
                if skymap is None:
                    skymap = self.getSkymap(namespace, "deepCoadd")

                for ref in namespace.butler.subset("calexp", dataId=dataId):
                    if not ref.datasetExists("calexp"):
                        continue

                    # XXX fancier mechanism to select an individual exposure than just pulling out "visit"?
                    visit = ref.dataId["visit"]
                    if visit not in visitRefs:
                        visitRefs[visit] = list()
                    visitRefs[visit].append(ref)

                    md = ref.get("calexp_md", immediate=True)
                    wcs = lsst.afw.image.makeWcs(md)
                    box = lsst.afw.geom.Box2D(lsst.afw.geom.Point2D(0, 0),
                                              lsst.afw.geom.Point2D(md.get("NAXIS1"), md.get("NAXIS2")))
                    # Going with just the nearest tract.  Since we're throwing all tracts for the visit
                    # together, this shouldn't be a problem unless the tracts are much smaller than a CCD.
                    tract = skymap.findTract(wcs.pixelToSky(box.getCenter()))
                    if overlapsTract(tract, wcs, box):
                        if visit not in visitTract:
                            visitTract[visit] = set()
                        visitTract[visit].add(tract.getId())
            else:
                tract = dataId.pop("tract")
                # making a DataRef for src fills out any missing keys and allows us to iterate
                for ref in namespace.butler.subset("src", dataId=dataId):
                    self._addDataRef(namespace, ref.dataId, tract)

        # Ensure all components of a visit are kept together by putting them all in the same set of tracts
        for visit, tractSet in visitTract.iteritems():
            for ref in visitRefs[visit]:
                for tract in tractSet:
                    self._addDataRef(namespace, ref.dataId, tract)
        if visitTract:
            tractCounter = collections.Counter()
            for tractSet in visitTract.itervalues():
                tractCounter.update(tractSet)
            log.info("Number of visits for each tract: %s" % (dict(tractCounter),))


def overlapsTract(tract, imageWcs, imageBox):
    """Return whether the image (specified by Wcs and bounding box) overlaps the tract

    @param tract: TractInfo specifying a tract
    @param imageWcs: Wcs for image
    @param imageBox: Bounding box for image
    @return bool
    """
    tractWcs = tract.getWcs()
    tractCorners = [tractWcs.pixelToSky(lsst.afw.geom.Point2D(coord)).getVector() for
                    coord in tract.getBBox().getCorners()]
    tractPoly = convexHull(tractCorners)

    try:
        imageCorners = [imageWcs.pixelToSky(lsst.afw.geom.Point2D(pix)) for pix in imageBox.getCorners()]
    except lsst.pex.exceptions.LsstCppException, e:
        # Protecting ourselves from awful Wcs solutions in input images
        if (not isinstance(e.message, lsst.pex.exceptions.DomainErrorException) and
            not isinstance(e.message, lsst.pex.exceptions.RuntimeErrorException)):
            raise
        return False

    imagePoly = convexHull([coord.getVector() for coord in imageCorners])
    if imagePoly is None:
        return False
    return tractPoly.intersects(imagePoly) # "intersects" also covers "contains" or "is contained by"


class ForcedPhotCcdConfig(ProcessImageForcedConfig):
    doApplyUberCal = lsst.pex.config.Field(
        dtype = bool,
        doc = "Apply meas_mosaic ubercal results to input calexps?",
        default = False
    )

## @addtogroup LSST_task_documentation
## @{
## @page processForcedCcdTask
## ForcedPhotCcdTask
## @copybrief ForcedPhotCcdTask
## @}

class ForcedPhotCcdTask(ProcessImageForcedTask):
    """!
    A command-line driver for performing forced measurement on CCD images

    This task is a subclass of ForcedPhotImageTask which is specifically for doing forced
    measurement on a single CCD exposure, using as a reference catalog the detections which
    were made on overlapping coadds.

    The run method (inherited from ForcedPhotImageTask) takes a lsst.daf.persistence.ButlerDataRef
    argument that corresponds to a single CCD.  This should contain the data ID keys that correspond to
    the "forced_src" dataset (the output dataset for ForcedPhotCcdTask), which are typically all those
    used to specify the "calexp" dataset (e.g. visit, raft, sensor for LSST data) as well as a coadd
    tract.  The tract is used to look up the appropriate coadd measurement catalogs to use as references
    (e.g. deepCoadd_src; see CoaddSrcReferencesTask for more information). While the tract must be given
    as part of the dataRef, the patches are determined automatically from the bounding box and WCS of the
    calexp to be measured, and the filter used to fetch references is set via config
    (BaseReferencesConfig.filter).

    In addition to the run method, ForcedPhotCcdTask overrides several methods of ForcedPhotImageTask
    to specialize it for single-CCD processing, including makeIdFactory(), fetchReferences(), and
    getExposure().  None of these should be called directly by the user, though it may be useful
    to override them further in subclasses.
    """

    ConfigClass = ForcedPhotCcdConfig
    RunnerClass = lsst.pipe.base.ButlerInitializedTaskRunner
    _DefaultName = "forcedPhotCcd"
    dataPrefix = ""

    def makeIdFactory(self, dataRef):
        """Create an object that generates globally unique source IDs from per-CCD IDs and the CCD ID.

        @param dataRef       Data reference from butler.  The "ccdExposureId_bits" and "ccdExposureId"
                             datasets are accessed.  The data ID must have the keys that correspond
                             to ccdExposureId, which is generally the same that correspond to "calexp"
                             (e.g. visit, raft, sensor for LSST data).
        """
        expBits = dataRef.get("ccdExposureId_bits")
        expId = long(dataRef.get("ccdExposureId"))
        return lsst.afw.table.IdFactory.makeSource(expId, 64 - expBits)

    def fetchReferences(self, dataRef, exposure):
        """Return a SourceCatalog of sources which overlap the exposure.

        The returned catalog is sorted by ID and guarantees that all included children have their
        parent included and that all Footprints are valid.

        @param dataRef       Data reference from butler corresponding to the image to be measured;
                             should have tract, patch, and filter keys.
        @param exposure      lsst.afw.image.Exposure to be measured (used only to obtain a Wcs and
                             bounding box).

        All work is delegated to the references subtask; see CoaddSrcReferencesTask for information
        about the default behavior.
        """
        references = lsst.afw.table.SourceCatalog(self.references.schema)
        badParents = set()
        unfiltered = self.references.fetchInBox(dataRef, exposure.getBBox(), exposure.getWcs())
        for record in unfiltered:
            if record.getFootprint() is None or record.getFootprint().getArea() == 0:
                if record.getParent() != 0:
                    self.log.warn("Skipping reference %s (child of %s) with bad Footprint" %
                    (record.getId(), record.getParent()))
                else:
                    self.log.warn("Skipping reference parent %s with bad Footprint" % (record.getId(),))
                    badParents.add(record.getId())
            elif record.getParent() not in badParents:
                references.append(record)
        references.sort()   # need to ensure catalog is in ID order so find methods work
        return references

    def getExposure(self, dataRef):
        """Read input exposure to measure

        @param dataRef       Data reference from butler.  Only the 'calexp' dataset is used,
                             unless config.doApplyUberCal is true, in which case the corresponding
                             meas_mosaic outputs are used as well.
        """
        exposure = ProcessImageForcedTask.getExposure(self, dataRef)
        if not self.config.doApplyUberCal:
            return exposure
        if applyMosaicResults is None:
            raise RuntimeError(
                "Cannot use improved calibrations for %s because meas_mosaic could not be imported."
                % dataRef.dataId
                )
        else:
            applyMosaicResults(dataRef, calexp=exposure)
        return exposure

    def _getConfigName(self):
        """!Return the name of the config dataset.  Forces config comparison from run-to-run
        """
        return self.dataPrefix + "forcedPhotCcd_config"

    def _getMetadataName(self):
        """!Return the name of the metadata dataset.  Forced metadata to be saved
        """
        return self.dataPrefix + "forcedPhotCcd_metadata"

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=PerTractCcdDataIdContainer)
        return parser


