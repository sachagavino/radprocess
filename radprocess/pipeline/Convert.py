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
from radprocess.constants.constants import Ggram, kbol, pc, amu


class Convert:

    def write_amr_to_zarr(output, store_path="amr_grid.zarr"):
        """
        args: the output of RAMSES converted by PyMses. 
        output: store in Zarr format the RAMSES outputs.
        """

        root = zarr.open(store_path, mode="w")

        N = output["x"].shape[0]

        # Scalars
        scalar_fields = [
            "Tgas", "Tdust", "x", "y", "z", "dx", "level", "density"
        ]

        for name in scalar_fields:
            root.create_dataset(
                name,
                data=output[name],
                chunks=(100_000,),
                compressor=zarr.Blosc(cname="zstd", clevel=3),
            )

        # Vector field
        root.create_dataset(
            "velocity",
            data=output["velocity"],
            chunks=(50_000, 3),
            compressor=zarr.Blosc(cname="zstd", clevel=3),
        )

        # Dust ratios
        root.create_dataset(
            "dust_ratio",
            data=output["dust_ratio"],
            chunks=(50_000, output["dust_ratio"].shape[1]),
            compressor=zarr.Blosc(cname="zstd", clevel=3),
        )

        root.attrs["n_cells"] = N

        return root






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
        has_velocity = any("velocity" in s for s in source)
        has_bfield = any("B_" in s for s in source)
        has_pressure = any("pressure" in s for s in source)
        has_temp = any("temperature" in s for s in source)

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
                dustratios[i-1,:] = output["dust_enrich{:d}".format(i)]/correction_factor
            
            #useful because ramses dumps the total density = rho_g + rho_d in units of mu micron.(so density is not the gas density)
            #for the dust, ramses actually provides dust enrichment (not dust-to-gas mass ratio)
 
            #create new array for dust-to-gas mass ratio

            # # correct dust to gas mass ratios
            # for i in range(1,nb_sizes+1):
            #     output[f'dust_enrich{(i + 1):d}'] = output["dust_enrich{:d}".format(i)]/correction_factor
            
        else:
            output["gas_massdensity"]  = output["density"]*mp # in unit of RAMSES here
            output['dust_massdensity'] = output['gas_massdensity'] * 1e-02
            

        #---VELOCITY SECTION:---
        if has_velocity:
            output['velocity'] = cells['velocity'] * snap.info['unit_velocity'].express(C.m / C.s)

        #---TEMPERATURE SECTION:---
        if has_temp:
            print("*************************************************")
            print("---------------ADDING TEMPERATURE----------------")
            # ============================================
            # Everything in cgs for the moment:
            # if 'mu_gas' in snap.info:
            #     mu = snap.info['mu_gas'] # mean molecular weight in amu
            # else:
            #     mu = 1.4
            ## mp       = mu * 1.660531e-24    # n gramme
            print("WORKING WITH mu_gas = ", mp/amu)
            scale_n  = 1.
            scale_l  = pc
            scale_d  = scale_n * mp
            scale_t  = 1.0 / np.sqrt(Ggram * scale_d)
            scale_v  = scale_l / scale_t    
            scale_T2 = mp/kbol * scale_v**2
            # ============================================
            unit = snap.info['unit_temperature'].express(C.K)
            print("unit", unit)
            print("scale_T2", scale_T2)
            # X = 0.76 # Hydrogen mass fraction
            output['Tgas'] = cells['thermal_pressure'] / cells['density'] * scale_T2
            output['Tdust'] = cells['thermal_pressure'] / cells['density'] * scale_T2
            print("Min Max Mean T", output['Tgas'].min(), output['Tgas'].max(), output['Tgas'].mean())
            print("*************************************************")

        # # correct dust to gas mass ratios
        # for i in range(1,nb_sizes+1):
        #     dustratios[i-1,:] = output["dustratio{:d}".format(i)]/correction_factor

        # # case for a single dust-to-gas mass ratio:
        # dustratio = np.sum(dustratios[:,:], axis=0)

        # print("output['Tgas']: ", output['Tgas'].shape)
        # print("output['Tdust']: ", output['Tdust'].shape)
        # print("output['x']: ", output['x'].shape)
        # print("output['y']: ", output['y'].shape)
        # print("output['z']: ", output['z'].shape)
        # print("output['dx']: ", output['dx'].shape)
        # print("output['level']: ", output['level'].shape)
        # print("output['density']: ", output['density'].shape)
        # print("output['velocity']: ", output['velocity'].shape)
        # print("output['dust_ratio']: ", output['dust_ratio'].shape)

        z_dustratios = zarr.array(dustratios)
        return z_dustratios # so it is stored in Grid, then used to write RT files.
        

    def polaris_photon(self):
        self.rat1 = 1

    def toramdc(self):
        self.rat2 = 1

    def topolaris(self):
        self.rat3 = 1

    