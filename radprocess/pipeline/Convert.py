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
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc

import pymses
from pymses.utils import constants as C
from pymses.filters import CellsToPoints

from radprocess import radmc3d
from radprocess import ramses

from radprocess.pipeline.OcTree import OcTree
from radprocess.constants.constants import Ggram, kbol, pc, amu


class Convert:
    def __init__(self):
        self.l_cm = None
        self.unit_l = None

    def write_amr_to_zarr(self,
                          output, 
                          store_path="amr_grid.zarr", 
                          ramses_num=None,
                          has_temp=False, 
                          has_velocity=False, 
                          has_ratio=False,
                          has_fluids=False,
                          has_fluid_v=False):
        """
        args: the output list of RAMSES converted by PyMses. 
        output: RAMSES outputs cleanly stored in Zarr format for efficient use.
        """

        root = zarr.open(store_path, mode="w")
        compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)

        N = output["x"].shape[0]

        # --Scalars--
        if has_temp:
            scalar_fields = [
                "Tgas", "Tdust", "x", "y", "z", "dx", "level", "gas_massdensity"
            ]
        else:
            scalar_fields = [
                "x", "y", "z", "dx", "level", "gas_massdensity"
            ]        

        for name in scalar_fields:
            arr = output[name]
            root.create_dataset(
                name,
                shape=arr.shape,
                dtype=arr.dtype,
                data=arr,
                chunks=(100_000,),
                compressor=compressor,
            )

            # root.create_dataset(
            #     name,
            #     data=output[name],
            #     chunks=(100_000,),
            #     compressor=Blosc(cname="zstd", clevel=3),
            # )

        # --Dust_massdensity-- (even though it is scalar we define a specific chunk because of possible multi-fluid model)
        arr = output["dust_massdensity"]

        if arr.ndim == 1:
            chunks = (50_000,)
        else:
            chunks = (50_000, arr.shape[1])

        root.create_dataset(
            "dust_massdensity",
            shape=arr.shape,
            dtype=arr.dtype,
            data=arr,
            chunks=chunks,
            compressor=compressor,
        )

        # --Velocity field--
        if has_velocity:
            arr = output["velocity"]
            root.create_dataset(
                "velocity",
                shape=arr.shape,
                dtype=arr.dtype,
                data=arr,
                chunks=(50_000, 3),
                compressor=compressor,
            )

        # --Fluid velocity field--
        if has_fluid_v:
            arr = output["fluid_v"]
            root.create_dataset(
                "fluid_velocity",
                shape=arr.shape,
                dtype=arr.dtype,
                data=arr,
                chunks=(50_000, 3),
                compressor=compressor,
            )

        if has_fluids:
            print('add here the grain sizes, and if exist, the dust number densities')

            # root.create_dataset(
            #     "dust_ratio",
            #     data=output["dust_ratio"],
            #     chunks=(50_000, output["dust_ratio"].shape[1]),
            #     compressor=Blosc(cname="zstd", clevel=3),
            # )


        # --Dust ratios--
        if has_ratio:
            arr = output["dust_ratio"]

            root.create_dataset(
                "dust_ratio",
                shape=arr.shape,
                dtype=arr.dtype,
                data=arr,
                chunks=(50_000, arr.shape[1]),
                compressor=compressor,
            )

        #root.attrs["n_cells"] = N

        root.attrs.update({
            "schema_version": int(1),
            "nb_cells": int(N),
            "l_cm": float(self.l_cm),
            "l_m": float(self.unit_l),
        })

        if ramses_num is not None:
            root.attrs["ramses_output_num"] = int(ramses_num)

        return root



    def ramses(self, ramses_dir, ramses_num, ramses_out, source, sim_param, nb_sizes):
        importlib.reload(pymses) # Reload pymses to clear internal caches. Another problem with Pymses.

        print(f"\n=== New RAMSES Conversion ===")
        print(f"RAMSES Folder: {ramses_dir}")
        print(f"Output: {ramses_num}")
        print(f"Dust species detected: {nb_sizes}")
        print(f"Enabled AMR fields: {source}\n")
    
        CLR_LINE = " " * 50 + "\r"
        cell_counter = 0
        fields = []
        i = 0
        #----------------------------------------------------
        


        #separate hydro vs. fluids
        hydro_fields = [s for s in source if not s.startswith("fluid_")]
        has_ratio    = any(s.startswith("dust_ratio") for s in source)
        has_velocity = any(s.startswith("velocity") for s in source)
        has_bfield   = any(s.startswith("B_") for s in source)
        has_pressure = any("pressure" in s for s in source)
        has_temp     = any("temperature" in s for s in source)


        fluid_fields = [s for s in source if s.startswith("fluid_")]
        has_fluids   = any(s.startswith("fluid_density") for s in source)
        has_fluid_v  = any(s.startswith("fluid_v") for s in source)

        snap = pymses.RamsesOutput(ramses_dir, ramses_num)

        # if nb_sizes > 0 and nb_sizes != snap.info["ndust"]:
        #     raise ValueError(
        #         "Mismatch between grain sizes in hydro_file_descriptor.txt "
        #         "and snap.info['ndust']"
        #     )

        snap.amr_fields()
        mp = snap.info["unit_density"].express(C.g_cc) #nu*mu in gcc

        # # Safest way to get the list:
        # available_fields = list(snap.info)
        # print("Pymses sees fields:", available_fields)

        # # Validate user-requested fields
        # for f in source:
        #     if f not in available_fields:
        #         raise RuntimeError(f"Requested field '{f}' not found in RAMSES data")

        
        ## HYDRO PART
        #amr = snap.amr_source(source)
        amr = snap.amr_source(hydro_fields)
        cell_source = CellsToPoints(amr)
        cell_source.ndim
        cells = cell_source.flatten()
        output = {} #OR = list size of amr + 4 for dx x y z level
        self.unit_l = snap.info["unit_length"].express(C.m) # Cell lengths 
        self.l_cm = snap.info["unit_length"].express(C.cm) # Cell lengths 
        # max. number of cells
        output["dx"] = cells.get_sizes()*self.unit_l
        nr_of_cells = len(output["dx"])
        # Original cell positions (from 0 to 1) converted into uint length)
        output["x"] = cells.points[:,0]*self.unit_l
        output["y"] = cells.points[:,1]*self.unit_l
        output["z"] = cells.points[:,2]*self.unit_l
        # level of each cell
        output["level"]=np.log2(self.unit_l/output["dx"])


        ## FLUID PART
        if fluid_fields:
            fluid_amr = snap.amr_source(fluid_fields)
            fluid_cells = CellsToPoints(fluid_amr).flatten()

        
        #

        # if has_ratio:
        #     dustratios = np.zeros((nb_sizes, nr_of_cells))
        #     for i in range(1,nb_sizes+1):
        #         output["dust_enrich{:d}".format(i)]  = cells["dust_ratio_{:d}".format(i)]  

        #     epsilon_tot = np.zeros(output["dust_enrich1"].shape)
        #     for i in range(1,nb_sizes+1):
        #         epsilon_tot+=output["dust_enrich{:d}".format(i)]
        #     correction_factor = 1-epsilon_tot
        #     output["gas_massdensity"] = mp*correction_factor*output["density"] # actual gas density 
        
        #     for i in range(1,nb_sizes+1):
        #         #dustratios[i-1,:] = output["dust_enrich{:d}".format(i)]/correction_factor
        #         output["dust_ratio{:d}".format(i)] = output["dust_enrich{:d}".format(i)]/correction_factor

        #---DENSITY SECTION:---
        #output["density"] = cells["density"]


        if has_ratio:
            #IMPORTANT: Here, ramses dumps the total density = rho_g + rho_d in units of mu micron.(so density is not the actual gas density)
            #for the dust, RAMSES provides dust enrichment (not dust-to-gas mass ratio). So we derive a corretion factor to obtain the dust mass ratios.

            # shape: (cells, species)
            dust_ratio = np.zeros((nr_of_cells, nb_sizes), dtype=np.float32)
            dust_massdensity = np.zeros((nr_of_cells, nb_sizes), dtype=np.float32)

            # RAMSES gives enrichment; convert to dust/gas
            epsilon_tot = np.zeros(nr_of_cells, dtype=np.float32)

            enrich = []

            for i in range(1, nb_sizes + 1):
                e = cells[f"dust_ratio_{i}"]
                enrich.append(e)
                epsilon_tot += e

            correction_factor = 1.0 - epsilon_tot
            output["gas_massdensity"] = mp * correction_factor * cells["density"]

            for i, e in enumerate(enrich):
                dust_ratio[:, i-1] = e / correction_factor


            for i in range(0, nb_sizes):
                dust_massdensity[:, i] = dust_ratio[:, i]*output["gas_massdensity"]

            output['dust_massdensities'] = dust_massdensity
            output['dust_massdensity']  = np.sum(dust_massdensity, axis=1)
            
        if has_fluids:
            fluid_density = np.zeros((nr_of_cells, nb_sizes), dtype=np.float32)

            for i in range(0, nb_sizes):
                fluid_density[:, i] = fluid_cells[f"fluid_density_{i+1}"] * mp

            output["dust_massdensities"] = fluid_density
            output['dust_massdensity']  = np.sum(fluid_density, axis=1)
            output["gas_massdensity"]  = cells["density"]*mp
            

        if not has_ratio and not has_fluids:
            output["gas_massdensity"]  = cells["density"] * mp # in unit of RAMSES here
            output['dust_massdensity'] = cells["density"] * mp * sim_param.dtogas
        
 
        #---VELOCITY SECTION:---
        if has_velocity:
            output["velocity"] = (
                cells["velocity"] * snap.info["unit_velocity"].express(C.m / C.s)
            )

        # #---FLUID VELOCITY SECTION:---
        # if has_fluid_v:
        #     output['fluid_v'] = fluid_cells['fluid_v']

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


        root = self.write_amr_to_zarr(output, 
                                 store_path = ramses_out / f"output{ramses_num}_grid.zarr", 
                                 ramses_num=ramses_num, 
                                 has_temp=has_temp, 
                                 has_velocity=has_velocity, 
                                 has_ratio=has_ratio,
                                 has_fluids=has_fluids,
                                 has_fluid_v=has_fluid_v)


        # print("shape of output['Tgas']: ", output['Tgas'].shape)
        # print("shape of output['x']: ", output['x'].shape)
        # print("shape of output['gas_massdensity']: ", output['gas_massdensity'].shape)
        # print("shape of output['dust_massdensity']: ", output['dust_massdensity'].shape)
        # print("shape of output['dust_massdensities']: ", output['dust_massdensities'].shape)
        # print("shape of output['fluid_v']: ", output['fluid_v'].shape)
        # print("shape of output['velocity']: ", output['velocity'].shape)
        # print(output['dust_massdensities'][0,0])




        return root # will later be stored in the main Grid, then used to write RT files.
        

    def to_radmc(self, root):

        l_m = root.attrs.get("l_m")
        nb_cells = root.attrs.get("nb_cells")
        level = root["level"]

        x_min = -0.5*l_m
        y_min = -0.5*l_m
        z_min = -0.5*l_m

        max_level = max(level)
        min_level = min(level)

        print("\n")
        print("Octree parameter:")
        print("    Level        (min,max)  : ", int(min_level),",", int(max_level))
        print("    Nr. of cells (data, max): ", nb_cells,",", 8**max_level)
        print("    Length       (min,max)  : ", l_m/(2**max_level),",", l_m, "\n")

        # Initialize the octree
        tree = OcTree.OcTree(x_min, y_min, z_min, l_m)


    def to_polaris(self):
        self.rat3 = 1


    def polaris_photon(self):
        self.rat1 = 1



    