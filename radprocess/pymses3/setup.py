#!/usr/bin/env python
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

import numpy

from Cython.Build import cythonize
from setuptools import setup, Extension, find_packages


def setup_pymses():
    setup(name="PyMSES",
          python_requires='>=3.6, <4',  # 3.6/3.7/3.8/3.9/3.10
          packages=find_packages(include=['pymses', 'pymses.*']),

          # Metadata to display on PyPI
          version="4.1.5",
          description="Analysis and visualization Python modules for RAMSES",
          # long_description="""TODO""",
          classifiers=["Intended Audience :: Science/Research",
                       "License :: OSI Approved :: GNU General Public License (GPL)",
                       "Operating System :: MacOS :: MacOS X",
                       "Operating System :: POSIX :: AIX",
                       "Operating System :: POSIX :: Linux",
                       "Programming Language :: C",
                       "Programming Language :: Python",
                       "Topic :: Scientific/Engineering :: Astrophysics",
                       "Topic :: Scientific/Engineering :: Analysis",
                       "Topic :: Scientific/Engineering :: Visualization",
                       "Topic :: Scientific/Engineering :: Post-processing"],
          keywords='astrophysics visualization amr ramses',
          author="Damien CHAPON, Thomas GUILLET, Olivier IFFRIG, Marc LABADENS",
          author_email="damien.chapon@cea.fr",
          url="http://irfu.cea.fr/Projets/PYMSES/index.html",
          license="GPL-3",
          setup_requires=['setuptools>=18.0',
                          'cython>=0.29.21',
          ],
          include_dirs=[numpy.get_include()],
          install_requires=["numpy>=1.19.1", ],
          ext_modules=[Extension('pymses.sources.ramses.tree_utils',
                                 sources=['pymses/sources/ramses/tree_utils.pyx']),
                       Extension('pymses.sources.ramses._octree_utils',
                                 sources=['pymses/sources/ramses/py_octree_utils.c']),
                       Extension('pymses.utils.point_utils',
                                 sources=['pymses/utils/point_utils.pyx']),
                       Extension('pymses.utils._point_utils',
                                 sources=['pymses/utils/py_point_utils.c']),
                       Extension('pymses.utils._ray_cube_utils',
                                 sources=['pymses/utils/py_ray_cube_utils.c']),
                       Extension('pymses.sources.ramses._read_ramses',
                                 sources=["pymses/sources/ramses/io/fio.c",
                                          "pymses/sources/ramses/io/read_amr.c",
                                          "pymses/sources/ramses/io/read_cells.c",
                                          "pymses/sources/ramses/io/read_parts.c",
                                          "pymses/sources/ramses/io/py_read_ramses.c"]),
                       Extension('pymses.sources.ramses.hilbert',
                                 sources=['pymses/sources/ramses/hilbert.pyx']),
                       Extension('pymses.analysis.raytracing._raytrace',
                                 sources=['pymses/analysis/raytracing/py_raytrace.c']),
                       Extension("pymses.analysis.ppv._ray_cast_ppv",
                                 sources=["pymses/analysis/ppv/py_ray_cast_ppv.c"]),
                       Extension("pymses.analysis.lic._line_integral_convolution",
                                 sources=["pymses/analysis/lic/py_line_integral_convolution.c"]),
          ],

    )
    return


if __name__ == '__main__':
    setup_pymses()
