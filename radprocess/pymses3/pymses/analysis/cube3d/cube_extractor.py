# -*- coding: utf-8 -*-
#   This file is part of PyMSES.
#
#   PyMSES is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PyMSES is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.
r"""
:mod:`pymses.analysis.cube3d.cube_extractor` --- 3D datacube extracction module
-------------------------------------------------------------------------------

"""
import numpy

from .cube import Datacube3D
from ..operator import ScalarOperator
from ..data_processor import DataProcessor, DataProcess


def amr2cube_process_task(idomain, amr_source, operator, points_datamap, points, interpolation, verbose):
    # Select points from current domain only
    ipoint_domain = numpy.nonzero(points_datamap == idomain)[0]
    points_in_domain = points[ipoint_domain, :]

    # Get tha AMR dataset #idomain
    amr_dset = amr_source.get_domain_dset(idomain, verbose)

    # Sample AMR fields in the domain #idomain
    points_dset = amr_dset.sample_points(points_in_domain, add_level=False, add_cell_center=False,
                                         interpolation=interpolation)

    return ipoint_domain, points_dset


class AmrToCubeProcess(DataProcess):
    def __init__(self, queues, amr_source, operator, points_datamap, points, interpolation, verbose):
        tasks_queue, results_queue = queues
        super(AmrToCubeProcess, self).__init__(tasks_queue, results_queue)
        self.amr_source = amr_source
        self.operator = operator
        self.points_datamap = points_datamap
        self.points = points
        self.interpolation = interpolation
        self.ipoint_domain_list = []
        self.points_dset_list = []
        self.verbose = verbose

    def process_task(self, idomain):
        task_res = amr2cube_process_task(idomain, self.amr_source, self.operator, self.points_datamap, self.points,
                                         self.interpolation, self.verbose)

        self.ipoint_domain_list.append(task_res[0])
        self.points_dset_list.append(task_res[1])
        return None

    @property
    def result(self):
        res = (self.ipoint_domain_list, self.points_dset_list)
        return res


class CubeExtractor(DataProcessor):
    r"""
    3D datacube extraction data processor class

    Parameters
    ----------
    source: : class:`~pymses.sources.ramses.sources.RamsesAmrSource`
        AMR data source.
    op : :class:`~pymses.analysis.operator.AbstractOperator`
            physical scalar quantity data operator
    verbose: ``bool``
        verbosity boolean flag. Default None.
    """

    def __init__(self, source, op, verbose=None):
        if not isinstance(op, ScalarOperator):
            raise NotImplementedError("3D datacube extraction is only available for a single scalar quantity so far.")
        super(CubeExtractor, self).__init__(source, op, amr_mandatory=True, verbose=verbose)
        self._interpolation = None
        self._points_datamap = None
        self._points = None

    def process_factory(self, tasks_queue, results_queue, verbosity):
        return AmrToCubeProcess((tasks_queue, results_queue), self._source, self._operator, self._points_datamap,
                                self._points, self._interpolation, verbosity)

    def process(self, camera, cube_size, resolution=256, interpolation=True):
        r"""
        3D cube extraction method : convert AMR octree datasets to cartesian regular grid.

        Parameters
        ----------
        camera : :class:`~pymses.analysis.camera.Camera`
            camera containing all the view params
        cube_size: ``float``
            3D datacube size (must be positive)
        resolution: `` int``
            datacube number of voxels along each dimension. Default 256.
        interpolation : ``boolean``
            CIC interpolation flage. Default True.

        Returns
        -------
        TODO
        """
        self._interpolation = interpolation
        self._points, required_level = camera.get_datacube_info(cube_size, resolution)

        # AMR DataSource preparation
        self._source.set_read_levelmax(required_level)

        # Compute the domain ids of the points
        if self._source.dom_decomp is not None:
            self._points_datamap = self._source.dom_decomp.map_points(self._points)
        else:
            self._points_datamap = numpy.ones(self._points.shape[0])
        # Sort the points by datamap
        unique_domains = numpy.unique(self._points_datamap)
        ndomains = unique_domains.size

        ipoint_batches = []
        point_dsets = []

        if self.use_multiprocessing and ndomains >= 2:
            self.init_multiprocessing_run(ndomains)

            # Loop over unique domains : enqueue jobs
            for idomain in unique_domains:
                self.enqueue_task(idomain)

            for res in self.get_results():
                ipoint_batches += res[0]
                point_dsets += res[1]
        else:
            # Loop over unique domains
            for idomain in unique_domains:
                task_res = amr2cube_process_task(idomain, self._source, self._operator, self._points_datamap,
                                                 self._points, self._interpolation, self._verbose)

                ipoint_batches.append(task_res[0])
                point_dsets.append(task_res[1])

        # Perform final reordering concatenation
        cat_method = point_dsets[0].concatenate
        dset = cat_method(point_dsets, reorder_indices=ipoint_batches)
        del point_dsets
        del ipoint_batches

        vdict = {}
        for key, func in self._operator.iter_scalar_func():
            vdict[key] = func(dset).reshape(resolution, resolution, resolution)
        del dset
        # Return single datacube with unit defined in the operator
        cube = Datacube3D(camera, cube_size, resolution, self._operator.output_unit)
        cube.data[...] = self._operator.operation(vdict)
        del vdict

        return cube


__all__ = ["CubeExtractor"]
