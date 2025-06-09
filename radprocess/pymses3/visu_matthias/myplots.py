import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import colormaps as cmaps
import quantity as Q
from functools import partial
from pymses.analysis import DataMap
from pymses.utils import constants as C
from pymses.utils.regions import Sphere
from pymses.analysis import (Camera, slicing, ScalarOperator,
                             FractionOperator, raytracing, bin_spherical)
from pymses.analysis.point_sampling import PointSamplingProcessor
from pymses.filters import (CellsToPoints, PointRandomDecimatedFilter,
                            PointFunctionFilter, RegionFilter)

GRAV_COMPAT=True

""" Contains several plot classes : IMF, histogram, hexbin, slice and columndensity"""

class myPlot(object):
    def __init__(self, output):
        self.output = output

    def _save(self, fig, fileName, fileFormat=".png"):
        savedir = os.path.join(self.output.pathName, "results")
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        fig.savefig(os.path.join(savedir, fileName+fileFormat))

class IMF(myPlot):
    plot_name = "IMF"

    def __call__(self, bins=10, range=None, save=False, ax=None, **kwargs):
        print("-> Entering %s ..." % self.plot_name)
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        ms = np.log10(self.output.sinkDataFrame["mass"].as_matrix())
        dn, bin_edge = np.histogram(ms, bins=bins, range=range)
        dndlogm = dn/np.diff(bin_edge)

        ax.step(bin_edge[1:], np.log10(dndlogm), where='pre', lw=1)
        ax.set(xlabel="log(mass)", ylabel="log(dn/d(logm))")

        if save:
            fileName = self._getFileName()
            self._save(fig, fileName, ".png")
            plt.close(fig)
        else:
            fig.show()

        return fig, ax

    def _getFileName(self):
        fileName = self.output.folderName+"_"+self.output.outputName
        fileName += "_"+self.plot_name
        return fileName

class hist(myPlot):
    plot_name = "hist"

    def __call__(self, x, bins=10, range=None, normed=False, weights=None,
                 cumulative=False, bottom=None, histtype='bar', align='mid',
                 orientation='vertical', rwidth=None, color=None, label=None,
                 stacked=False, logx=False, logy=False, ax=None, save=False, **kwargs):
        print("-> Entering %s ..." % self.plot_name)
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        # Loads needed fields in the _qtyDataFrame
        if weights is None:
            fields_to_unload = self.output.load(x)
            ylabel = "Cells"
            wweights = None
        else:
            fields_to_unload = self.output.load(x, weights)
            ylabel = weights.label()
            wweights = self.output.qtyDataFrame[weights.name].as_matrix()
        xx = self.output.qtyDataFrame[x.name].as_matrix()
        # Unloads fields that have been loaded
        self.output.unload(*fields_to_unload)

        # Handles log or linear scale
        if logx:
            xscale = "log"
            if isinstance(bins, int):
                if range is None:
                    bins = np.logspace(np.log10(xx.min()), np.log10(xx.max()),
                                       bins)
                else:
                    bins = np.logspace(np.log10(range[0]), np.log10(range[1]),
                                       bins)
        else:
            xscale = "linear"

        # Plots histogram
        ax.hist(xx, bins, range, normed, wweights, cumulative, bottom,
                histtype, align, orientation, rwidth, logy, color,
                label, stacked, **kwargs)
        ax.set(xlabel=x.label(), ylabel=ylabel, xscale=xscale)

        # Defines title
        info = self.output.info
        title = "%.1d kyr" % (info["time"] * info["unit_time"].express(C.kyr))
        if self.output.massAcc is not None:
            title += " (mass accreted : %3.1d %%)" % (self.output.massAcc*100)
        ax.set_title(title)
        fig.tight_layout()

        if save:
            fileName = self._getFileName(x, weights)
            self._save(fig, fileName, ".png")
            plt.close(fig)
        else:
            fig.show()

        return fig, ax

    def _getFileName(self, x, weights):
        fileName = self.output.folderName+"_"+self.output.outputName
        fileName += "_"+self.plot_name+"_"+x.name
        if weights is not None:
            fileName += "_weights_"+weights.name
        return fileName

class hexbin(myPlot):
    plot_name = "hexbin"

    def __call__(self, x, y, c=None, gridsize=300, bins=None, xscale='lin',
                 yscale='lin', extent=None, cmap=cmaps.viridis,
                 norm=colors.LogNorm(), vmin=None, vmax=None, alpha=None,
                 edgecolors='none', reduce_C_function=np.sum, mincnt=1e-300,
                 ax=None, save=False, **kwargs):
        print("-> Entering %s ..." % self.plot_name)
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        # Loads needed fields in the _qtyDataFrame
        if c is None:
            fields_to_unload = self.output.load(x, y)
            cc = None
            cblabel = "Cells"
        else:
            fields_to_unload = self.output.load(x, y, c)
            cc = self.output.qtyDataFrame[c.name].as_matrix()
            cblabel = c.label()
        xx = self.output.qtyDataFrame[x.name].as_matrix()
        yy = self.output.qtyDataFrame[y.name].as_matrix()

        # Unloads fields that have been loaded
        self.output.unload(*fields_to_unload)

        # Plots hexbin
        im = ax.hexbin(xx, yy, cc, gridsize, bins, xscale, yscale, extent,
                       cmap, norm, vmin, vmax, alpha, edgecolors=edgecolors,
                       reduce_C_function=reduce_C_function, mincnt=mincnt, **kwargs)
        cb = ax.get_figure().colorbar(mappable=im, label=cblabel)

        ax.set(xlabel=x.label(), ylabel=y.label())

        # Defines title
        info = self.output.info
        title = "%.1d kyr" % (info["time"] * info["unit_time"].express(C.kyr))
        if self.output.massAcc is not None:
            title += " (mass accreted : %3.1d %%)" % (self.output.massAcc*100)
        ax.set_title(title)

        if save:
            fileName = self._getFileName(x, y, c)
            self._save(fig, fileName, ".png")
            plt.close(fig)
        else:
            fig.show()

        return fig, ax, cb

    def _getFileName(self, x, y, c):
        fileName = self.output.folderName+"_"+self.output.outputName
        fileName += "_"+self.plot_name+"_"+x.name+"_"+y.name
        if c is not None:
            fileName += "weights_"+weigths.name
        return fileName

class MeanProfile(myPlot):
    plot_name = "MeanProfile"
    def __call__(self, typlot, radius, ax=None, save=False, **kwargs):
        print("-> Entering %s ..." % self.plot_name)
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        liste = []
        dxmin = self.output.cellSize("max",None)

        for Id in self.output.sinkDataFrame["Id"]:
            print(Id)
            center = self.output.getCenter(Id,None)
            r_bins, profile = self.smooth_log_profile(self.output.ramsesOutput,
                                                       typlot, center,
                                                       radius, dxmin, **kwargs)
            bins_center = (r_bins[:-1]+r_bins[1:])/2
            ax.plot(bins_center, profile, 'k', lw=1)
            liste.append(profile)
        print("")

        listprofile = np.array(liste)
        minprofile = np.min(listprofile, axis=0)
        maxprofile = np.max(listprofile, axis=0)
        meanprofile = np.average(listprofile,
                                 weights=self.output.sinkDataFrame["mass"],
                                 axis=0)
        varprofile = np.average((meanprofile-listprofile)**2,
                                weights=self.output.sinkDataFrame["mass"],
                                axis=0)
        stdprofile = np.sqrt(varprofile)
        
        ax.plot(bins_center, minprofile, 'k', lw=1)
        ax.plot(bins_center, maxprofile, 'k', lw=1)
        ax.fill_between(bins_center, meanprofile-stdprofile, meanprofile+stdprofile, alpha=0.4)
        ax.plot(bins_center, meanprofile, 'b', lw=2)
        ax.set(xscale='log', yscale='log')
        ax.set(xlabel=r"$\rm r \ \left( \ UA \ \right)$", ylabel=typlot.label())

        # Defines title
        title = "%.1d kyr" % (self.output.info["time"] * 
                              self.output.info["unit_time"].express(C.kyr))
        if self.output.massAcc is not None:
            title += " (mass accreted : %3.1d %%)" % (self.output.massAcc*100)
        ax.set_title(title)
        if save:
            fileName = self._getFileName(typlot)
            self._save(fig, fileName, ".png")
            plt.close(fig)
        else:
            fig.show()
        return fig, ax

    @staticmethod
    def _sph_profile(ro, typlot, center, radius, bins, verbose=False, npoints=1.0e7):
        source = ro.amr_source(typlot.amrFields,
                               grav_compat=GRAV_COMPAT,
                               verbose=False)
        sphere = Sphere(center,radius)
        points = sphere.random_points(npoints)
        point_dset = PointSamplingProcessor(source).process(points, add_level=True)
        del points

        weight_func = typlot.get_calc(ro)
        profile = bin_spherical(point_dset, center, weight_func,
                                bins, divide_by_counts=True)

        return profile

    @classmethod
    def smooth_log_profile(cls, ro, typlot, center, rMax, rMin, **kwargs):
        r_bins = np.logspace(np.log10(rMin),np.log10(rMax),100)
        profile = cls._sph_profile(ro, typlot, center, rMax, r_bins, **kwargs)

        radius_c = 8*rMin
        if rMax > radius_c:
            mask = r_bins<radius_c
            r_bins2 = r_bins[mask]
            profile2 = cls._sph_profile(ro, typlot, center, 2*radius_c, r_bins2, npoints=1e7)
            profile[:len(profile2)] = profile2
        r_bins *= ro.info["unit_length"].express(C.au)
        return r_bins, profile

    def _getFileName(self, typlot):
        fileName = self.output.folderName+"_"+self.output.outputName
        fileName += "_"+self.plot_name+"_"+typlot.name
        return fileName


##########################################################################
# Map2d plots
##########################################################################

class map2d(myPlot):
    """
    Abstrac class
    Common interface to slice and column plots
    """
    def __call__(self, typlot, vector_field=None, 
                 center=[0.5, 0.5, 0.5], los="x", size=[0.5, 0.5], up_vector=None,
                 log_sensitive=True, distance=0.5, far_cut_depth=0.5,
                 map_max_size=512, add_sink_info=False, save=False,
                 axis_unit=C.pc, **kwargs):
        print("-> Entering %s ..." % self.plot_name)
        output = self.output
        cam = Camera(output.getCenter(center, None), los, up_vector, size,
                     output.info["unit_length"],
                     distance, far_cut_depth, map_max_size, log_sensitive)
        
        if add_sink_info:
            sinkdset=output.sinkDataFrame[["x","y","z","mass"]]
        else:
            sinkdset=None

        # Creates datamap
        mymapx = self.getDataMap(ro=output.ramsesOutput, typlot=typlot,
                                 vector_field=vector_field,
                                 sinkdset=sinkdset,
                                 cam=cam, verbose=output._verbose)
        fig = mymapx.save_plot(axis_unit=axis_unit, **kwargs)

        # Saving
        if save:
            fname = self._getFileName(typlot, center, size, los)
            self._save(fig, plot_fname, ".png")
            plt.close(fig)
        else:
            fig.show()
        return fig

    def _getFileName(self, typlot, center, size, los):
        if isinstance(center, str):
            if center == "sink_barycenter":
                centerName = "C_sink_barycenter"
        elif isinstance(center, int):
            centerName = "C_sinkId%d" % center
        else:
            centerName = "C_%02d_%02d_%02d" % (10*center[0],10*center[1],10*center[2])
        sizeName = "S_%02d_%02d" % (10*size[0], 10*size[1])
        if isinstance(los, str):
            losName = "L_%s" % los
        else:
            losName = "L_%02d_%02d_%02d" % (10*los[0],10*los[1],10*los[2])
        fileName = self.output.folderName+"_"+self.output.outputName
        fileName += "_"+self.plot_name+"_"+typlot.name
        fileName += "_"+losName+"_"+centerName+"_"+sizeName
        return fileName

class slice(map2d):
    plot_name = "slice"

    # Creates a pymses DataMap and then converts it into a personal datamap to handle sink info
    # and vector fields
    @staticmethod
    def getDataMap(ro, typlot, vector_field=None, sinkdset=None, cam=None, verbose=False):
        """ sinkdset is an array containing all the data needed to plot sink info"""
        if cam is None:
            cam = Camera(size_unit=ro.info["unit_length"])
        # Gets list of amr fields needed for amr_source of pymses
        source = ro.amr_source(typlot.amrFields,
                               grav_compat=GRAV_COMPAT,
                               verbose=verbose)

        # Defines Operator of pymses needed for the slicing
        varplot_op = ScalarOperator(typlot.get_calc(ro), typlot.unit)
        
        # Slicing process
        mapx = slicing.SliceMap(source, cam, varplot_op, z=0.0)

        if vector_field is not None:
            source = ro.amr_source(vector_field.amrFields,
                                   grav_compat=GRAV_COMPAT,
                                   verbose=verbose)
            mapx._vector_field = vector_field
            u_axis, v_axis, _ = cam.get_camera_axis()
            
            vect_U = vector_field.dot(u_axis)
            varplot_op = ScalarOperator(vect_U.get_calc(ro), vect_U.unit)
            mapxU = slicing.SliceMap(source, cam, varplot_op, z=0.0)
            mapx._vmapU = mapxU._vmap.copy().transpose()

            vect_V = vector_field.dot(v_axis)
            varplot_op = ScalarOperator(vect_V.get_calc(ro), vect_V.unit)
            mapxV = slicing.SliceMap(source, cam, varplot_op, z=0.0)
            mapx._vmapV = mapxV._vmap.copy().transpose()
        else:
            mapx._vector_field = None
        
        mapx._cblabel = typlot.label()
        mapx._ro = ro
        mapx._info = ro.info
        mapx._sinkdset = sinkdset
        mapx = myDataMap.convert2myDataMap(mapx)
        return mapx

class column(map2d):
    plot_name = "column"

    # Creates a pymses DataMap and then converts it into a personal datamap to handle sink info
    # and vector fields
    @staticmethod
    def getDataMap(ro, typlot, vector_field=None, sinkdset=None, cam=None, verbose=False):
        if cam is None:
            cam = Camera(size_unit=ro.info["unit_length"])

        # Gets list of amr fields needed for amr_source of pymses
        amr_fields_to_read = Q.get_amrfields_to_read(typlot, Q.rho)
        source = ro.amr_source(amr_fields_to_read,
                               grav_compat=GRAV_COMPAT,
                               verbose=verbose)
        if typlot.name == Q.rho.name:
            func = partial(lambda dset, func, factor: func(dset)*factor,
                           func=Q.rho.get_calc(ro, unit=C.g/C.cm**3),
                           factor=ro.info['unit_length'].express(C.cm))
            unit = Q.rho.unit * ro.info['unit_length']
            varplot_op = ScalarOperator(func, unit)
            cblabel = r"$n \ ( \ g.cm^{-2} \ )$"
        else:
            up_func = partial(lambda dset, f, g: f(dset)*g(dset),
                              f=Q.rho.get_calc(ro),
                              g=typlot.get_calc(ro))
            down_func = Q.rho.get_calc(ro)
            frac_unit = typlot.unit
            varplot_op = FractionOperator(up_func, down_func, frac_unit)
            cblabel = typlot.label()

        rt = raytracing.RayTracer(source, ro.info, varplot_op)
        mapx = rt.process(cam, surf_qty=True)

        if vector_field is not None:
            amr_fields_to_read = Q.get_amrfields_to_read(vector_field, Q.rho)
            source = ro.amr_source(amr_fields_to_read,
                                   grav_compat=GRAV_COMPAT,
                                   verbose=verbose)
            u_axis, v_axis, _ = cam.get_camera_axis()
            down_func = Q.rho.get_calc(ro)
            frac_unit = vector_field.unit
            
            vect_U = vector_field.dot(u_axis)
            up_func = partial(lambda dset, f, g: f(dset)*g(dset),
                              f=Q.rho.get_calc(ro),
                              g=vect_U.get_calc(ro))
            varplot_op = FractionOperator(up_func, down_func, frac_unit)
            rt = raytracing.RayTracer(source, ro.info, varplot_op)
            mapxU = rt.process(cam, surf_qty=True)

            vect_V = vector_field.dot(v_axis)
            up_func = partial(lambda dset, f, g: f(dset)*g(dset),
                              f=Q.rho.get_calc(ro),
                              g=vect_V.get_calc(ro))
            varplot_op = FractionOperator(up_func, down_func, frac_unit)
            rt = raytracing.RayTracer(source, ro.info, varplot_op)
            mapxV = rt.process(cam, surf_qty=True)

            mapx._vmapU = mapxU._vmap.copy().transpose()
            mapx._vector_field = vector_field
            mapx._vmapV = mapxV._vmap.copy().transpose()
        else:
            mapx._vector_field = None

        # Adds some info needed to DataMap object
        mapx._cblabel = cblabel
        mapx._ro = ro
        mapx._info = ro.info
        mapx._sinkdset = sinkdset
        mapx = myDataMap.convert2myDataMap(mapx)

        return mapx

########################################################################
# Subclass of the pymses class "DataMap"
########################################################################
class myDataMap(DataMap):
    @classmethod
    def convert2myDataMap(cls, datamap):
        datamap.__class__ = myDataMap
        return datamap

    # Redefinition of the save_plot function of the DataMap class
    def save_plot(self, vrange=None, fraction=1.0, cmap="Viridis",
                  discrete=False, is_log_values=None, axis_unit=None,
                  map_unit=None, verbose=False, **kwargs):
        fig = super(myDataMap, self).save_plot(None, vrange, fraction, cmap,
                                               discrete, is_log_values, axis_unit,
                                               map_unit, verbose)
        ax, cbax = fig.get_axes()

        # Sets ylabel for the colorbar
        cbax.set_ylabel(self._cblabel)

        # Adds vector field
        if self._vector_field is not None:
            self._add_vector_field(axis_unit, ax, **kwargs)
        
        # Adds sink informations on the plot
        if self._sinkdset is not None:
            arraySinkPos = np.array([self._sinkdset["x"],
                                     self._sinkdset["y"],
                                     self._sinkdset["z"]])
            arraySinkPos = arraySinkPos.T
            arraySinkMass = np.array(self._sinkdset["mass"])
            self._add_sink_info(arraySinkPos, arraySinkMass, axis_unit, ax)
        
        # Adds title
        time = self._info["time"] * self._info["unit_time"].express(C.kyr)
        title = "%.3f kyr" % time
        ax.set_title(title)

        fig.tight_layout()

        return fig

    # Adds sink info on the "ax" axis
    def _add_sink_info(self, arraySinkPos, arraySinkMass, axis_unit, ax):
        # Gets edges of each pixel of the camera in axis_unit
        uvaxes_edges_labels = self.camera.get_uvaxes_edges_labels(axis_unit)
        (u_axisname, _, _, uedges), (v_axisname, _, _, vedges) = uvaxes_edges_labels
        # Gets vectors of the orthonormal base defining the camera plane
        u_axis, v_axis, z_axis = self.camera.get_camera_axis()

        # Defines the accretion radius
        dxmin = 1. / 2**self._info["levelmax"]
        if axis_unit is not None:
            dxmin *= self._info["unit_length"].express(axis_unit)
        if("ncell_racc" in self._info):
            r_acc = self._info['ncell_racc'] * dxmin
        else:
            r_acc = 4 * dxmin

        if uedges[-1] < uedges[0]:
            u_axis = -u_axis
        if vedges[-1] < vedges[0]:
            v_axis = -v_axis

        sinkMassTot = np.sum(arraySinkMass)
        for (sinkPos, sinkMass) in zip(arraySinkPos, arraySinkMass):
            if u_axisname == 'u':
                up = np.dot(u_axis, sinkPos/self._info['boxlen'] - self.camera.center)
            else:
                up = np.dot(u_axis, sinkPos/self._info['boxlen'])
            if v_axisname == 'v':
                vp = np.dot(v_axis, sinkPos/self._info['boxlen'] - self.camera.center)
            else:
                vp = np.dot(v_axis, sinkPos/self._info['boxlen'])
            if axis_unit is not None:
                up *= self._info["unit_length"].express(axis_unit)
                vp *= self._info["unit_length"].express(axis_unit)
            col = cmaps.viridis(sinkMass/sinkMassTot)
            ax.add_patch(plt.Circle((up, vp), radius=r_acc, fc='white', alpha=0.3))
            ax.add_patch(plt.Circle((up, vp), radius=dxmin, fc=col, ec=col))

        ax.text(0.02, -0.09, "$\mathrm{M}_*=$ %4.1f $M_\odot$" % sinkMassTot,
                transform=ax.transAxes, color='k')

    # Adds vector field on the "ax" axis
    def _add_vector_field(self, axis_unit, ax, arrows=32, scale="max"):
        # Gets edges of each pixel of the camera in axis_unit
        uvaxes_edges_labels = self.camera.get_uvaxes_edges_labels(axis_unit)
        (_, _, _, uedges), (_, _, _, vedges) = uvaxes_edges_labels

        # Computes centers of pixels in axis_unit
        ucenters = (uedges[:-1] + uedges[1:])/2
        vcenters = (vedges[:-1] + vedges[1:])/2

        # Samples arrows according arrows parameter
        step = len(ucenters) / arrows
        start = step / 2 - 1
        U = self._vmapU[start::step, start::step]
        V = self._vmapV[start::step, start::step]
        ucenters = ucenters[start::step]
        vcenters = vcenters[start::step]

        # Defines grid to plot the sample of arrows
        X, Y = np.meshgrid(ucenters, vcenters)

        # Defines scaling for ploting arrows
        if scale == "max":
            scale = np.max(np.sqrt(U**2+V**2))
        elif scale == "mean":
            scale = np.sqrt(np.mean(U**2+V**2))

        Quiver = ax.quiver(X, Y, U, V, units='inches',scale=scale*3,
                           pivot='mid', color='w')
        ax.quiverkey(Quiver, 0.70, -0.08, scale,
                     r'$%.2e \ %s$' % (scale, self._vector_field.texUnitName),
                     labelpos='E', coordinates='axes', color='w', labelcolor='k')
        ax.axis([uedges[0], uedges[-1], vedges[0], vedges[-1]])
