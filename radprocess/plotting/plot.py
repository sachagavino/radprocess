#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
file name: plot
@author: Sacha Gavino
last update: June 2025
language: > PYTHON 3.8
__________________________________________________________________________________________
short description:  plotting routines
__________________________________________________________________________________________
"""
import glob, sys

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from ..constants.constants import au2cm, M_sun, R_sun



def density2D(mass1, mass2=None, overlap=False):
    grid = np.loadtxt('thermal/amr_grid.inp', engine='python', skiprows=5)
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    dens = np.loadtxt('thermal/dust_density.inp', engine='python', header=None, skiprows=3)
    dens = dens[0].values
    nbspecies = int(len(dens)/(nr*nt))
    dens = np.reshape(dens, (nbspecies, nt, nr))
    dist = grid[:nr+1]/autocm
    theta = grid[nr+1:nr+1+nt+1]
    theta[-1] = np.pi
    dist, tt = np.meshgrid(dist, theta)
    rr = dist*np.sin(tt)
    zz = dist*np.cos(tt)

    dens[dens<=1e-100] = 1e-100
    if overlap == False:
        # #--PLOT FIGURE--
        if nbspecies == 1:
            fig = plt.figure(figsize=(10, 8.))
            ax = fig.add_subplot(111)
            plt.xlabel(r'r [au]', fontsize = 17)
            plt.ylabel(r'z [au]', fontsize = 17, labelpad=-7.4)
            
            numdens = dens[0]#/mass1[0]
            t = plt.pcolor(rr, zz, numdens, cmap='gnuplot2', shading='auto', norm=LogNorm(vmin=1e-80, vmax=1e-1))
            clr = plt.colorbar(t)
            clr.set_label(r'$n_\mathrm{d}$ [cm${-3}$]', labelpad=-33, y=1.06, rotation=0, fontsize = 16)
            #plt.xlim(1, 500)
            #plt.ylim(-300, 300)
            ax.tick_params(labelsize=17)
            clr.ax.tick_params(labelsize=16) 
            plt.show()

        else:
            for ispec in range(nbspecies):
                fig = plt.figure(figsize=(8, 8.))
                ax = fig.add_subplot(111)
                plt.xlabel(r'r [au]', fontsize = 17)
                plt.ylabel(r'z [au]', fontsize = 17, labelpad=-7.4)
                numdens = dens[ispec]#/mass1[ispec]
                t = plt.pcolor(rr, zz, numdens, cmap='gnuplot2', shading='auto', norm=LogNorm(vmin=1e-30, vmax=1e-17))
                clr = plt.colorbar(t)
                clr.set_label(r'$\rho_\mathrm{d}$ [g.cm${-3}$]', labelpad=-33, y=1.06, rotation=0, fontsize = 16)
                plt.xlim(1, 1000)
                plt.ylim(-1000, 1000)
                ax.tick_params(labelsize=17)
                clr.ax.tick_params(labelsize=16) 
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
                ax.text(0.90, 0.95, 'bin: {}'.format(ispec+1), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=16, bbox=props)
                plt.show()

    if overlap == True:
        density = np.zeros((nt, nr))
        fig = plt.figure(figsize=(10, 10.))
        ax = fig.add_subplot(111)
        plt.xlabel(r'r [au]', fontsize = 17)
        plt.ylabel(r'z [au]', fontsize = 17, labelpad=-7.4)
        numdens = dens#/mass1[ispec]
        for ispec in range(0, nbspecies):
            density += numdens[ispec]
        t = plt.pcolor(rr, zz, density, cmap='gnuplot2', shading='auto', norm=LogNorm(vmin=1e-25, vmax=1e-17))
        clr = plt.colorbar(t)
        clr.set_label(r'$\rho_\mathrm{d}$ [g.cm${-3}$]', labelpad=-33, y=1.06, rotation=0, fontsize = 16)
        #plt.xlim(1, 500)
        #plt.ylim(-300, 300)
        ax.tick_params(labelsize=17)
        clr.ax.tick_params(labelsize=16) 
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.90, 0.95, 'env+disk', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=16, bbox=props)
        plt.show()
