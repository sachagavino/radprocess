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
:mod:`pymses.analysis.slicing.processor` --- 2D slice map processing module
---------------------------------------------------------------------------

"""
from ..point_sampling import PointSamplingProcessor
from ..data_processor import DataProcessor
from ..datamap import DataMap
from ..operator import MultiFieldOperator
from pymses.utils.constants import none as None_unit


class SliceMapProcessor(DataProcessor):
    r"""
    2D slice map processor class. Compute 2D slice maps out of an AMR data source

    Parameters
    ----------
    amr_source   : :class:`~pymses.sources.ramses.output.RamsesAmrSource`
        AMR data source
    operator : :class:`~pymses.analysis.operator.AbstractOperator`
        Data slicing operator.
    """
    def __init__(self, amr_source, operator):
        super(SliceMapProcessor, self).__init__(amr_source, operator, amr_mandatory=True)

    def process(self, camera, z=0.0, interpolation=False):
        """
        Compute a slice map made of sampling points.

        Parameters
        ----------
        camera : :class:`~pymses.analysis.camera.Camera`
            camera handling the view parameters
        z      : ``float``
            position of the slice plane along the line-of-sight axis of the camera
        interpolation : ``boolean``
            Need to perform interpolation in amr tree point sampling ? default False
            THIS IS NOT IMPLEMENTED YET.

        Returns
        -------
        map : ``array``
            sliced map
        """
        # Map size
        nx, ny = camera.get_map_size()  # Sampling points
        p = camera.get_slice_points(z)

        # Get dataset
        if self._source.ndim == 2:
            p = p[:, 0:2]  # 3D -> 2D restriction

        add_level = self._operator.use_cell_dx()
        max_lev = camera.get_required_resolution()

        # Point sampling processing
        psp = PointSamplingProcessor(self._source)
        if not self._user_needs_multiprocessing:
            psp.disable_multiprocessing()
        dset = psp.process(p, add_level=add_level, max_search_level=max_lev, interpolation=interpolation)

        # Compute map values according to operator
        maps = {}
        for key, func in self._operator:
            if add_level:
                # for MaxLevelOperator
                maps[key] = dset["level"].reshape(nx, ny)
            else:
                maps[key] = func(dset).reshape(nx, ny)

        if isinstance(self._operator, MultiFieldOperator):
            # Build datamap dictionary with consistent units
            map_dict = self._operator.operation(maps)
            output_unit_dict = self._operator.output_unit

            datamap = None
            for key, map in map_dict.iter_items():
                if key not in output_unit_dict:
                    print("Warning : undefined Unit for map key '%s'. Set to 'none' Unit." % key)
                    u = None_unit
                else:
                    u = output_unit_dict[key]
                if datamap is None:
                    datamap = DataMap(map, camera, u, map_name=key)
                else:
                    datamap.add_scalar_map(map, key, u)
        else: # Return single datamap with unit defined in the operator
            map = self._operator.operation(maps)
            datamap = DataMap(map, camera, self._operator.output_unit)

        return datamap


def SliceMap(source, camera, op, z=0.0, interpolation=False, use_multiprocessing=True):
    """
    Shortcut method for the :class:`~pymses.analysis.slicing.SliceMapProcessor` processing method.

    See Also
    --------
    processor : :class:`~pymses.analysis.slicing.SliceMapProcessor` class.
    proc_method : :meth:`~pymses.analysis.slicing.SliceMapProcessor.process` method.
    """
    p = SliceMapProcessor(source, op)
    if not use_multiprocessing:
        p.disable_multiprocessing()
    return p.process(camera, z, interpolation)


__all__ = ["SliceMapProcessor", "SliceMap"]
