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
"""
tree_utils.pyx -- low-level tree utility functions
"""

import numpy
cimport numpy
cimport cython

cdef inline double fmax (double x, double y):
    if x>y:return x
    return y
cdef inline double fmin (double x, double y):
    if x<y:return x
    return y
cdef extern from "math.h":
    double fabs(double x)
    double floor(double x)
    double ceil(double x)

@cython.wraparound(False)
@cython.boundscheck(False)
def tree_search(amr_struct, numpy.ndarray[double, ndim=2]points, max_search_level=None):
    """ Searches deepest cells containing points in a RAMSES tree structure.

    This function searches the deepest cells containing the points `points`,
    and returns the corresponding grid and cell indices.

    Parameters:

        amr_struct -- a dictionary tuple containing the AMR structure, such as
            the output of read_ramses_amr_file()[1]

        points -- a (npoints, ndim) numpy array holding the coordinates of the
            points to search for, which should be in [0, 1[^ndim

        max_search_level -- the deepest level to explore (int or a (npoints,) numpy array of int

    Output:

        a dictionary containing info about the matching cells.
    """

    points = numpy.asarray(points)

    cdef int npoints = points.shape[0]
    cdef char ndim = points.shape[1]
    cdef int max_read_level = amr_struct["readlmax"]

    cdef char c_max_search_level, c_max_temp
    cdef numpy.ndarray[char, ndim=1] c_max_search_level_array
    if max_search_level is None:
        c_max_search_level_array = max_read_level * \
                numpy.ones(npoints, 'i1')
    else:
        if isinstance(max_search_level, numpy.ndarray):
            assert max_search_level.size == npoints
            c_max_search_level_array = max_search_level.astype('i1')
        if isinstance(max_search_level, list):
            assert len(max_search_level)==npoints
            c_max_search_level_array = numpy.asarray(max_search_level, 'i1')
        if isinstance(max_search_level, int):
            c_max_search_level_array = max_search_level * \
                numpy.ones(npoints, 'i1')


    # Big AMR arrays
    cdef numpy.ndarray[double, ndim=2] xgrids = amr_struct["grid_centers"]
    cdef numpy.ndarray[int, ndim=2] sons = amr_struct["son_indices"]

    # Output arrays
    cdef numpy.ndarray[int, ndim=1] point_igrids = numpy.empty(npoints, 'i')
    cdef numpy.ndarray[char, ndim=1] point_icells = numpy.empty(npoints, 'i1')
    cdef numpy.ndarray[char, ndim=1] point_ilevels = numpy.empty(npoints, 'i1')

    cdef int ipoint, igrid, igridstart
    cdef char ibit, idim, ind, ilevel

    # This was in io_ramses.f90, but does not work
    #igridstart = amr_struct["coarse_son_indices"][0]
    igridstart = 0

    # Loop over points
    for ipoint in range(npoints):
        c_max_temp = c_max_search_level_array[ipoint]
        if c_max_temp>max_read_level:
            c_max_search_level = max_read_level
        else:
            c_max_search_level = c_max_temp

        # Start search with the topmost grid

        # Was in ioramses.f90: igrid = region_sons[0]
        igrid = igridstart

        ilevel = 1
        # Search down levels
        for ilevel in range(1, c_max_search_level+1):
            # Find the subcell
            ind = 0
            for idim in range(ndim):
                ibit = (points[ipoint, idim] > xgrids[igrid, idim])*1
                ind += (ibit << idim)

            # Is the subcell a leaf cell?
            if (sons[igrid, ind] < 0) or (ilevel == c_max_search_level):
                # We're done searching: the subcell is the one we want
                break
            else:
                # Go down the tree
                igrid = sons[igrid, ind]

        point_igrids[ipoint] = igrid
        point_icells[ipoint] = ind
        point_ilevels[ipoint] = ilevel

    return { "grid_indices" : point_igrids,
             "cell_indices" : point_icells,
             "levels" : point_ilevels }


@cython.wraparound(False)
@cython.boundscheck(False)
def cell_centers_from_grid_centers(numpy.ndarray[double, ndim=2] grid_centers, numpy.ndarray[int, ndim=1] grid_levels):
    """Compute the (ngrid, twotondim, ndim) array of the coordinates of the
    centers of the RAMSES cells from the grid centers and corresponding grid
    levels
    """

    cdef int ngrids = grid_centers.shape[0]
    cdef int ndim = grid_centers.shape[1]
    cdef int twotondim = 1 << ndim

    # Allocate output array
    cdef numpy.ndarray[double, ndim=3] cell_centers = \
            numpy.empty([ngrids, twotondim, ndim])

    # Precompute the cell shifts
    cdef numpy.ndarray[double, ndim=2] xc = numpy.empty([twotondim, ndim])
    cdef int idim, icell, ind, igrid
    cdef int ishift
    for ind in range(twotondim):
        for idim in range(ndim):
            ishift = (ind >> idim) & 1
            xc[ind, idim] = ishift - 0.5

    # Compute cell centers
    cdef double dx
    for igrid in range(ngrids):
        for ind in range(twotondim):
            for idim in range(ndim):
                dx = 0.5**grid_levels[igrid]
                cell_centers[igrid, ind, idim] = grid_centers[igrid, idim] \
                        + dx * xc[ind, idim]

    return cell_centers


@cython.wraparound(False)
@cython.boundscheck(False)
def grid_levels(amr_struct):
    levels = numpy.empty(amr_struct["ngrids"], 'i')

    # At each level, we store the grids of all CPUs, then all boundaries
    ngrids_level = amr_struct["ngridlevel"].sum(axis=0) + \
            amr_struct["ngridbound"].sum(axis=0)

    offset = 0
    for ilevel in range(len(ngrids_level)):
        ngrids = ngrids_level[ilevel]
        if ngrids > 0:
            levels[offset:(offset+ngrids)] = ilevel + 1
            offset += ngrids

    return levels


@cython.boundscheck(False)
@cython.wraparound(False)
def octree_build(octree_dset, amr_dset, camera_info, radius, include_split_cells=True):
    """ octree_build function which fills an octree structure with a RAMSES AMR dataset
    by taking into account the original AMR refinement only in a specific region.
    """
    # Big RAMSES AMR arrays
    cdef numpy.ndarray[int, ndim=2] ngridlevel = amr_dset.amr_struct["ngridlevel"]
    cdef numpy.ndarray[int, ndim=2] ngridbound = amr_dset.amr_struct["ngridbound"]
    cdef numpy.ndarray[int, ndim=2] son_indices = amr_dset.amr_struct["son_indices"]
    cdef numpy.ndarray[double, ndim=2] grid_centers = amr_dset.amr_struct["grid_centers"]
    cdef int idata=0
    cdef int readlmax = amr_dset.amr_struct["readlmax"]
    cdef int ngrids = amr_dset.amr_struct["ngrids"]
    cdef char ndim = amr_dset.amr_header["ndim"]
    cdef int icpu = amr_dset.icpu
    cdef char twotondim = 1 << ndim
    cdef char include_splt_cells = include_split_cells

    # RAMSES AMR structure igrid offset
    cdef int offset = 0
    # AMR cell size
    cdef double dx = 1.0 # Cell size (ilevel = 0 => dx=1.0)
    cdef double dx_2

    cdef char ind, idim, idim2, ishift, last_level, grid_ok, is_radius=0
    cdef int nbefore, nkeep, igrid, igrid_min, igrid_max, ngridl, ilevel, level
    cdef numpy.ndarray[double, ndim=2] xshift = numpy.empty((twotondim, ndim))
    cdef numpy.ndarray[char, ndim=1] ind_ok = numpy.zeros(twotondim, 'i1')
    cdef numpy.ndarray[double, ndim=2] cam_axes
    cdef numpy.ndarray[double, ndim=1] cell_box_min, cell_box_max, cam_box_min, cam_box_max, cam_center
    cdef numpy.ndarray[double, ndim=1] cell_center = numpy.zeros(ndim)
    cdef numpy.ndarray[double, ndim=1] cell_center2 = numpy.zeros(ndim)
    # New octree vars
    cdef int o_ngrid_max, o_igrid_max, o_igrid, o_igrid2, o_ilevel, o_level
    cdef char o_ind, o_ind2, o_ibit, o_idim
    cdef double o_dx, o_xshift, rad2=0.0, cc, cc2, ccmin, ccmax
    cdef numpy.ndarray[int, ndim=2] o_son_indices = octree_dset.amr_struct["son_indices"]
    cdef numpy.ndarray[double, ndim=2] o_grid_centers = octree_dset.amr_struct["grid_centers"]
    cdef numpy.ndarray[int, ndim=1] o_levels = octree_dset.amr_struct["cell_levels"]

    # Grid and cell masks for data gathering
    cdef numpy.ndarray[int, ndim=1] igrid_filter      = numpy.zeros(ngrids*twotondim, 'i')
    cdef numpy.ndarray[int, ndim=1] o_igrid_filter    = numpy.zeros(ngrids*twotondim, 'i')
    cdef numpy.ndarray[char, ndim=1] indcell_filter   = numpy.zeros(ngrids*twotondim, 'i1')

    ###### Precompute the cell shifts #######
    for ind in range(twotondim):
        for idim in range(ndim):
            ishift = (ind >> idim) & 1
            xshift[ind, idim] = ishift - 0.5
    #########################################

    # Camera info
    cam_center, cam_axes, cam_map_box, cell_box = camera_info
    cam_box_min, cam_box_max = cam_map_box.get_bounding_box()
    cell_box_min, cell_box_max = cell_box.get_bounding_box()
    if radius is not None:
        is_radius=1
        rad2 = radius*radius

    o_ngrid_max = octree_dset.amr_struct["ngrids"]
    o_igrid_max = octree_dset.amr_struct["o_igrid_max"]
    ##################################  Browse RAMSES AMR octree - AMR level iteration  ###########################################
    offset = 0
    for ilevel in range(readlmax): # Iterates over the AMR levels
        # Grid total (all cpus) number at current level
        ngridl = ngridlevel[:,ilevel].sum()+ngridbound[:,ilevel].sum()
        if ngridl==0:
            continue
        # Min. and max. active grid index for the current level
        nbefore = ngridlevel[:icpu, ilevel].sum()
        nkeep = ngridlevel[icpu, ilevel]
        if nkeep==0:
            offset += ngridl
            continue
        igrid_min = offset + nbefore
        igrid_max = igrid_min + nkeep
        # Cell level and size
        level = ilevel + 1
        last_level = (level==readlmax)
        dx = dx / 2.0
        dx_2 = dx / 2.0
        for igrid in range(igrid_min, igrid_max):# Iterates over RAMSES active grids
            # Iterates over RAMSES active cells and checks whether it contains leaf cells
            # intersecting with the camera/region
            grid_ok=0
            for ind in range(twotondim):
                if not include_splt_cells and (not last_level and son_indices[igrid, ind]>=0):
                    ind_ok[ind] = 0 # Flag non leaf-cell
                else:
                    # If include_splt_cells == True, all active cells are processed
                    # as if they were leaf cells
                    ind_ok[ind] = 1 # Flag leaf cell
                    grid_ok=1 # If there is at least one leaf cell, flag grid

            if not grid_ok: continue # The current grid does not contains leaf cell(s)

            # Check whether the leaf cell(s) intersects the camera/region
            grid_ok = 0
            if is_radius: # Detect whether the cell volume intersects a spherical region centered
            # on the camera center and of radius 'rad'
                for ind in range(twotondim):
                    if not ind_ok[ind]: continue # Skip non-leaf cell
                    cc2=0.0
                    for idim in range(ndim):
                        # Coordinate of the cell center from the camera center
                        cc = grid_centers[igrid, idim] + dx * xshift[ind, idim] - cam_center[idim]
                        # Adjusted coordinate of the nearest corner of the cell from the camera center
                        ccmin = fmin(fabs(cc - dx_2), fabs(cc + dx_2))
                        # Cum. of the square distance of the nearest cell corner
                        cc2 += ccmin*ccmin
                    if cc2 <= rad2: grid_ok=1 # Nearest cell corner in the sphere => intersecting cell !
                    else: ind_ok[ind] = 0 # Flag non-intersecting cells
            else: # Detect whether the cell volume intersects the camera region
                for ind in range(twotondim):
                    if not ind_ok[ind]: continue # Skip non-leaf cell
                    for idim in range(ndim):# Coordinate of the cell center from the camera center
                        cell_center[idim] = grid_centers[igrid, idim] + dx * xshift[ind, idim] - cam_center[idim]
                    # Cell center coordinates in the (u, v, w) camera system
                    for idim in range(ndim):
                        cell_center2[idim] = 0.0
                        for idim2 in range(ndim):
                            cell_center2[idim] += cam_axes[idim, idim2] * cell_center[idim2]

                        # Cell faces min./max. coordinates in the (u,v,w) camera system
                        ccmin = cell_center2[idim] + dx * cell_box_min[idim]
                        ccmax = cell_center2[idim] + dx * cell_box_max[idim]
                        if((ccmin >= cam_box_max[idim]) or (ccmax <= cam_box_min[idim])): # Non-intersecting cell
                            ind_ok[ind]=0 # Flag non-instersecting cell
                            break
                    if ind_ok[ind]:
                        grid_ok=1 # If at least one leaf cell intersects the camera region, flag grid

            if not grid_ok: continue # The current grid does not contain leaf cell(s) that intersect the camera/sphere region.

            # Tree-search in octree to find where to insert the current RAMSES grid
            o_igrid = 0
            o_dx = 1.0
            for o_ilevel in range(ilevel):
                o_level = o_ilevel + 1
                o_dx = o_dx / 2.0
                # Find son index in octree
                o_ind = 0
                for o_idim in range(ndim):
                    o_ibit = (grid_centers[igrid, o_idim] > o_grid_centers[o_igrid, o_idim])
                    o_ind += (o_ibit << o_idim)

                o_igrid2 = o_son_indices[o_igrid, o_ind]
                if o_igrid2<0:# Non-refined grid => create twotondim son grids in octree
                    o_levels[o_igrid_max] = o_level+1
                    # Set new grid coordinates
                    for o_idim in range(ndim):
                        o_xshift = xshift[o_ind, o_idim] * o_dx
                        o_grid_centers[o_igrid_max, o_idim] = o_grid_centers[o_igrid, o_idim] + o_xshift
                    o_son_indices[o_igrid, o_ind] = o_igrid_max
                    o_igrid2 = o_igrid_max
                    # Incremented octree grid index
                    o_igrid_max += 1
                    if o_igrid_max == o_ngrid_max:
                        # Record new octree current grid index
                        octree_dset.amr_struct["o_igrid_max"] = o_igrid_max
                        return

                o_igrid = o_igrid2

            # Deal with son cells of the grid
            for o_ind in range(twotondim):# Iterates over intersecting RAMSES active leaf cells
                if not ind_ok[o_ind]: continue
                igrid_filter[idata]   = igrid
                o_igrid_filter[idata] = o_igrid
                indcell_filter[idata] = o_ind
                idata += 1

        # End of RAMSES active grids iteration
        offset += ngridl # Cum. grid offset for next level
    ############################################# End of AMR levels iteration #####################################################

    # Record new octree current grid index
    octree_dset.amr_struct["o_igrid_max"] = o_igrid_max

    # Gather data
    for field in amr_dset.scalars:
        octree_dset[field][o_igrid_filter, indcell_filter] = amr_dset[field][igrid_filter, indcell_filter]
    for field in amr_dset.vectors:
        octree_dset[field][o_igrid_filter, indcell_filter, :] = amr_dset[field][igrid_filter, indcell_filter, :]
    for field in amr_dset.multivalued:
        octree_dset[field][o_igrid_filter, indcell_filter, :] = amr_dset[field][igrid_filter, indcell_filter, :]

    return


@cython.boundscheck(False)
@cython.wraparound(False)
def octree_compute_neighbors(octree_dset, verbose=True):
    """
    Find neighbour grids in octree

    """
    cdef char ndim = octree_dset.amr_header["ndim"]
    cdef char twotondim = 1 << ndim
    cdef char twondim = 2*ndim
    cdef int nreadlmax = octree_dset.amr_struct["readlmax"]
    cdef int ngrids = octree_dset.amr_struct["ngrids"]
    cdef int igrid, ind_son, son, ilevel, nson
    cdef int ileveln, igridn, indn
    cdef char idim, iface, isgn, sgn_face

    cdef numpy.ndarray[int, ndim=2] sons = octree_dset.amr_struct["son_indices"]
    cdef numpy.ndarray[int, ndim=2] nbor = octree_dset.amr_struct["neighbors"]

    # Tree piles
    cdef numpy.ndarray[int, ndim=1] igrid_nav = numpy.zeros(nreadlmax, dtype="i")
    cdef numpy.ndarray[int, ndim=1] ind_nav = numpy.zeros(nreadlmax, dtype="i")

    if verbose: print("Computing neighbor grid indices in octree")
    ilevel=0
    igrid = igrid_nav[ilevel]
    ind_son = ind_nav[ilevel]
    while ((ilevel != 0) or (ind_son != twotondim)):
        son = sons[igrid, ind_son]
        if son>0 and son<ngrids: # Compute son grid neighbors
            ilevel += 1
            igrid_nav[ilevel] = son
            ind_nav[ilevel] = -1
            for idim in range(ndim):
                for sgn_face in range(-1,2,2):
                    isgn = (-sgn_face+1)>>1
                    ileveln = ilevel - 1
                    igridn = igrid
                    indn = ind_son
                    ############################### Find neighbour grid in octree #############################
                    while ( (((indn>>idim)&1) != isgn)  and (ileveln>0) ): # Search the octree up to find neighbouring grid
                        ileveln -= 1
                        igridn = igrid_nav[ileveln]
                        indn = ind_nav[ileveln]

                    # Boundary grid : no neighbor
                    if ( (ileveln==0) and (((indn>>idim)&1) != isgn) ):continue

                    # Neighbour cell index in the current grid
                    indn += sgn_face * (1<<idim) # neighbour cell index

                    # Octree search down to the leaf-cell
                    while (ileveln < ilevel):
                        ileveln += 1
                        # Fetch son grid
                        nson = sons[igridn, indn]
                        if nson<0:break
                        igridn = nson
                        # Neighbour cell index
                        indn = ind_nav[ileveln] - sgn_face * (1<<idim)

                    # Sanity check : neighbour cells should never be more than 1 level of
                    # refinement apart (RAMSES octree property).
                    #if (ileveln < (ilevel-1)):
                    #	print "level warning : (ilevel, ilevel_neighbour) = ", ilevel, ileveln
#					else:
                    if (ileveln >= (ilevel-1)):
                        # Set neigbor grid indices
                        iface = 1 + 2 * idim - isgn
                        nbor[son, iface] = igridn
                    ###########################################################################################

        ind_son = ind_nav[ilevel] + 1
        while ind_son == twotondim:
            if ilevel==0:break
            ilevel -= 1
            ind_son = ind_nav[ilevel] + 1
        ind_nav[ilevel] = ind_son
        igrid = igrid_nav[ilevel]

    if verbose: print("Computing neighbor grid indices in octree : done")
    return
