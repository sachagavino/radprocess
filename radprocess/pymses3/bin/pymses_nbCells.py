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

# Count the number of leaf_mask cells of a ramses 3D output
CELL_LEVEL_HISTOGRAM = True
from time import time
t0 = time()
from pymses.utils.misc import chunks
from pymses import RamsesOutput
from numpy import repeat,newaxis,arange
try: #Use MPI if possible:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    MPI_process_number = comm.Get_size()
    if myrank == 0:print(("Import time =", time() - t0, "MPI_process_number = ", MPI_process_number))
except:
    myrank = 0
    MPI_process_number = 1
from optparse import OptionParser
parser = OptionParser()
parser.usage = "%prog ramses_directory ramses_output_number"
(opts, args) = parser.parse_args()
try:
    fileDir = args[0]
    outNumber = int(args[1])
except:
    # "None" leads to an automatic look for
    # a RAMSES output in the current directory
    fileDir = None
    outNumber = None
ro = RamsesOutput(fileDir, outNumber)
outNumber = ro.iout
ncpu = ro.info["ncpu"]
if myrank == 0: print(("levelmin =", ro.info["levelmin"],"levelmax =", ro.info["levelmax"]))
cpu_full_list = list(range(1,ncpu+1))
cpu_node_list = chunks(cpu_full_list, MPI_process_number)
ro.verbose = False
rho = ro.amr_source([])
nbCells = 0
leaf_cells_level = {}
t1 = time()
if myrank < len(cpu_node_list):
    for icpu in cpu_node_list[myrank]:
        dset = rho.get_domain_dset(icpu)
        grid_mask = dset.get_active_mask()
        cell_lvl = repeat(dset.get_grid_levels(grid_mask)[:, newaxis], \
            8, axis=1).reshape((-1,))
        lvlmax = dset.amr_struct["readlmax"]
        sons = dset.amr_struct["son_indices"][grid_mask, :].reshape((-1,))
        leaf_mask = (sons < 0)+(cell_lvl == lvlmax)
        nbCells += leaf_mask.sum()
        if CELL_LEVEL_HISTOGRAM:
            for ilvl in range(lvlmax):
                leaf_mask = (sons < 0)*(cell_lvl == ilvl)
                if ilvl in leaf_cells_level:
                    leaf_cells_level[ilvl] += leaf_mask.sum()
                else :
                    leaf_cells_level[ilvl] = leaf_mask.sum()
            ilvl = lvlmax
            leaf_mask = cell_lvl == lvlmax
            if ilvl in leaf_cells_level:
                leaf_cells_level[ilvl] += leaf_mask.sum()
            else :
                leaf_cells_level[ilvl] = leaf_mask.sum()
    if myrank == 0:
        for i in range(1, MPI_process_number):
            nbCells_I = comm.recv(source=i, tag=3)
            nbCells += nbCells_I
            if CELL_LEVEL_HISTOGRAM:
                leaf_cells_level_I = comm.recv(source=i, tag=2)
                for ilvl in list(leaf_cells_level_I.keys()):
                    if ilvl in leaf_cells_level:
                        leaf_cells_level[ilvl] += leaf_cells_level_I[ilvl]
                    else:
                        leaf_cells_level[ilvl] = leaf_cells_level_I[ilvl]
        print(("Processing time =", time() - t1))
        print(("Total time =", time() - t0))
        print(("nbCells =",nbCells,"~ %.3e in output"%(nbCells), ro.output_repos,\
            ro.iout, "with ncpu =", ncpu))
        if CELL_LEVEL_HISTOGRAM:
            print(leaf_cells_level)
            pos = arange(len(list(leaf_cells_level.keys())))
            width = 1.0     # gives histogram aspect to the bar diagram
            import matplotlib.pyplot as plt
            ax = plt.axes()
            ax.set_xticks(pos + (width / 2))
            ax.set_xticklabels(list(leaf_cells_level.keys()))
            plt.bar(pos, list(leaf_cells_level.values()), width)
#			plt.show()
            plt.savefig("leaf_tree_cell_lvl_number_histogram_%s.png"%(outNumber))
            assert(sum(leaf_cells_level.values()) == nbCells)
    else:
        comm.send(nbCells, dest=0, tag=3)
        if CELL_LEVEL_HISTOGRAM:
            comm.send(leaf_cells_level, dest=0, tag=2)


