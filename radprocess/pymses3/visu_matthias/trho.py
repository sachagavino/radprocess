#!/usr/import/epd/python
# -*- coding: utf-8 -*-
#
# Script from M. Gonzalez
#
# run trho.py -d ./ -p ./ -o 89 
import sys
import getopt
import pymses as pm
import pylab  as pl
import numpy  as np
import os
from pymses.utils import constants as ct

bold  = "\033[1m"
reset = "\033[0;0m"

def usage():
    print((bold + "Usage: %%s [<options>]" + reset))
    print((bold + "-h --help:      " + reset + "print usage summary"))
    print((bold + "-p --path:      " + reset + "path of directory containing the outputs"))
    print((bold + "-d --directory: " + reset + "directory of the outputs"))
    print((bold + "-o --output:    " + reset + "number of the output to plot"))
    print((bold + "-r --record:    " + reset + "save the plot"))

def main(argv):
    try:
        opts,args = getopt.getopt(argv, "hp:d:o:",
           ["help", "path=", "directory=", "output=","record"])
    except getopt.GetoptError:
        usage()
        sys.exit(0)

    path      = "./"
    directory = "./"
    output    = 1
    record    = False

    for o, a in opts:
        if o in ["-h", "help"]:
            usage()
            sys.exit(0)
        elif o in ["-p", "path"]:
            path = a
        elif o in ["-d", "directory"]:
            directory = a
        elif o in ["-o", "output"]:
            output = eval(a)
        elif o in ["-r", "record"]:
            record = True

    myplot(path=path,directory=directory,output=output,record=record)


def myplot(path='./',directory='./',output=1,record=False):
    
    filerep       = path + directory
    if(record):
        savedir       = filerep+"results/trho/"
        if not os.path.exists(savedir):
            os.makedirs(savedir)
    
    ro = pm.RamsesOutput(filerep, output)

    time   = ro.info["time"]*ro.info["unit_time"].express(ct.kyr)
    
    # Define the output data structure (to use with pymses5)
    ro.define_amr_scalar_field("hydro", "rho", 0)
    ro.define_amr_vector_field("hydro", "vel", [1, 2, 3])
    #ro.define_amr_multivalued_field("hydro", "B", 4, 6)
    ro.define_amr_vector_field("hydro", "Bl", [4, 5, 6])
    ro.define_amr_vector_field("hydro", "Br", [7, 8, 9])
    ro.define_amr_scalar_field("hydro", "P", 10)
    ro.define_amr_multivalued_field("hydro", "Er", 11, ro.info["ngrp"])
    ro.define_amr_vector_field("hydro", "J", [11+ro.info["ngrp"], 12+ro.info["ngrp"], 13+ro.info["ngrp"]])
    ro.define_amr_scalar_field("hydro", "eint", 14+ro.info["ngrp"])
    if ro.info["eos"]:
        ro.define_amr_scalar_field("hydro", "T", 15+ro.info["ngrp"])
    ro.define_amr_vector_field("grav", "g", [0, 1, 2])

    # Read the fields
    if ro.info["eos"]:
        source = ro.amr_source(["rho","T","Bl","Br"])
    else:
        source = ro.amr_source(["rho","P","Bl","Br"])

    csource = pm.filters.CellsToPoints(source)
    cells   = csource.flatten()
    rho     = cells.fields['rho']*ro.info['unit_density'].express(ct.g/(ct.cm)**3)
    if ro.info["eos"]:
        T = cells.fields['T']
    else:
        T = cells.fields['P']/cells.fields['rho']*ro.info["unit_temperature"].express(ct.K)*ro.info["mu_gas"]
    # B = 1./4*( (cells.fields['Bl'][:,0]+cells.fields['Br'][:,0])**2  \
    #           +(cells.fields['Bl'][:,1]+cells.fields['Br'][:,1])**2  \
    #           +(cells.fields['Bl'][:,2]+cells.fields['Br'][:,2])**2)
    # B = B * 4*np.pi*(ro.info['unit_density']*ro.info['unit_velocity']**2).express(ct.erg/(ct.cm)**3)
    # print min(B),max(B)

    dx=cells.get_sizes()
    mass=cells.fields['rho']*ro.info['unit_density'].express(ct.Msun/(ct.cm)**3)*(dx*ro.info['unit_length'].express(ct.cm))**3

    fig = pl.figure()

    pl.xlabel(r'log($\rho$) [g/cm$^3$]')
    pl.ylabel(r'log(T) [K]')

    Tmin=0.95*min(np.log10(T))
    Tmax=1.05*max(np.log10(T))
    rhomin=1.05*min(np.log10(rho)) #because negative...
    rhomax=0.95*max(np.log10(rho)) #because negative...

    H,yedges,xedges=np.histogram2d(np.log10(T),np.log10(rho),bins=500,weights=mass,range=[[Tmin,Tmax],[rhomin,rhomax]])
    pl.imshow(np.log10(H),interpolation='none',origin='low',extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],aspect='auto')
    pl.colorbar()

    #pl.ylabel(r'Magnetic field (Gauss)')
    #pl.loglog(rho,B,'b.',ms=5.)
    #pl.axis([1e-19,1e-9,1.e-5,1.e1])
    #pl.scatter(rho,T,s=np.ones(np.size(rho)),c='r')

    pl.title("t=%.3f kyr" %time)

    if(record):
        for format in ("png", "eps", "pdf"):
            namefig = savedir + "trho" + "_%05d." %output + format
            print((bold + "Save figure: " + reset + namefig))
            pl.savefig(namefig)
        pl.close()
    else:
        pl.show(block=False)

if __name__ == "__main__":
    main(sys.argv[1:])
