"""
_____________________________________________________________________________________________________________
file name: Convert
last update: June 2025
language: > PYTHON 3.8
short description: convert from ramses to radmc3d and polaris. 
_____________________________________________________________________________________________________________
"""
import os
import importlib

import numpy as np
import zarr

import pymses
from pymses.utils import constants as C
from pymses.filters import CellsToPoints

from radprocess import radmc3d
from radprocess import ramses
from radprocess.pipeline.OcTree import OcTree
import radprocess.constants.constants as cst


class Convert:

    def ramses(self, ramses_folder, ramses_num, radmc_dir, source, sim_param, nb_sizes):
        importlib.reload(pymses) # Reload pymses to clear internal caches. Another problem with Pymses.

        print(f"\n=== New RAMSES Conversion ===")
        print(f"Folder: {ramses_folder}")
        print(f"Output: {ramses_num}")
        print(f"Dust species detected: {nb_sizes}")
        print(f"Enabled AMR fields: {source}\n")
    
        CLR_LINE = " " * 50 + "\r"
        cell_counter = 0
        fields = []
        i = 0
        #----------------------------------------------------
        has_dust = any("dust_ratio" in s for s in source)
        snap = pymses.RamsesOutput(ramses_folder, ramses_num)

        if nb_sizes > 0 and nb_sizes != snap.info["ndust"]:
            raise ValueError(
                "Mismatch between grain sizes in hydro_file_descriptor.txt "
                "and snap.info['ndust']"
            )

        snap.amr_fields()
        mp = snap.info["unit_density"].express(C.g_cc) #nu*mu in gcc

        # # Safest way to get the list:
        # available_fields = list(snap.info)
        # print("Pymses sees fields:", available_fields)

        # # Validate user-requested fields
        # for f in source:
        #     if f not in available_fields:
        #         raise RuntimeError(f"Requested field '{f}' not found in RAMSES data")

        amr = snap.amr_source(source)
        cell_source = CellsToPoints(amr)
        cell_source.ndim
        cells = cell_source.flatten()
        output = {} #OR = list size of amr + 4 for dx x y z level
        unit_l = snap.info["unit_length"].express(C.m) # Cell lengths 
        # max. number of cells
        output["dx"] = cells.get_sizes()*unit_l
        
        nr_of_cells = len(output["dx"])
        # Original cell positions (from 0 to 1) converted into uint length)
        output["x"] = cells.points[:,0]*unit_l
        output["y"] = cells.points[:,1]*unit_l
        output["z"] = cells.points[:,2]*unit_l
        # level of each cell
        output["level"]=np.log2(unit_l/output["dx"])

        #---DENSITY SECTION:---
        output["density"]  = cells["density"]
        if has_dust:
            dustratios = np.zeros((nb_sizes, nr_of_cells))
            for i in range(1,nb_sizes+1):
                output["dust_enrich{:d}".format(i)]  = cells["dust_ratio_{:d}".format(i)]  

            epsilon_tot = np.zeros(output["dust_enrich1"].shape)
            for i in range(1,nb_sizes+1):
                epsilon_tot+=output["dust_enrich{:d}".format(i)]
            correction_factor = 1-epsilon_tot
            output["gas_massdensity"] = mp*correction_factor*output["density"] # actual gas density 
        
            for i in range(1,nb_sizes+1):
                dustratios[i-1,:] = cells["dust_ratio_{:d}".format(i)]/correction_factor
            
            #useful because ramses dumps the total density = rho_g + rho_d in units of mu micron.(so density is not the gas density)
            #for the dust, ramses actually provides dust enrichment (not dust-to-gas mass ratio)
 
            #create new array for dust-to-gas mass ratio

            # # correct dust to gas mass ratios
            # for i in range(1,nb_sizes+1):
            #     output[f'dust_enrich{(i + 1):d}'] = output["dust_enrich{:d}".format(i)]/correction_factor
            
        else:
            output["gas_massdensity"]  = output["density"]*mp # in unit of RAMSES here
            output['dust_massdensity'] = output['gas_massdensity'] * 1e-02
            

        # #create new array for dust-to-gas mass ratio
        # dustratios = np.zeros((nb_sizes, nr_of_cells))

        # # correct dust to gas mass ratios
        # for i in range(1,nb_sizes+1):
        #     dustratios[i-1,:] = output["dustratio{:d}".format(i)]/correction_factor

        # # case for a single dust-to-gas mass ratio:
        # dustratio = np.sum(dustratios[:,:], axis=0)

        # #scales
        # scale_n = 1.
        # scale_l = cst.pc2cm
        # scale_d = scale_n*mp
        # scale_t = 1.0/np.sqrt(cst.Ggram*scale_d)
        # scale_v = scale_l / scale_t    
        # scale_T2 = mp/cst.kbol * scale_v**2

        #add condition regarding the content of source (T, P, velocity, etc.)


        z_dustratios = zarr.array(dustratios)
        return z_dustratios #radiative friendly format for inside Pipeline.py, so it is stored in Grid, then used to write radmc3d.inp files.
        

    def ramses2polaris(self):
        self.rat3 = 1

    