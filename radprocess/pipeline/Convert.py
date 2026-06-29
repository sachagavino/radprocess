"""
_____________________________________________________________________________________________________________
file name: Convert
last update: Feb 2026
language: > PYTHON 3.9
short description: convert from ramses to radmc3d and polaris. 
_____________________________________________________________________________________________________________
"""
import os
import sys
import re
from collections import defaultdict
import importlib
from pathlib import Path

import numpy as np
import zarr
from zarr.codecs import BloscCodec

import pymses
from pymses.utils import constants as C
from pymses.filters import CellsToPoints

from radprocess import radmc3d
from radprocess import ramses

from radprocess.pipeline.OcTree import OcTree, CellOct
from radprocess.constants.constants import Ggram, kbol, pc2cm, pc2m, amu, au2m, au2cm, M_sun, L_sun, R_sun, sigma


class Convert:
    def __init__(self):
        self.l_cm = None
        self.unit_l = None

    def write_amr_to_zarr(self,
                          output,
                          nb_species,
                          store_path="amr_grid.zarr", 
                          ramses_num=None,
                          has_temp=False, 
                          has_velocity=False, 
                          has_ratio=False,
                          has_fluids=False,
                          has_fluid_v=False,
                          has_bfield=False):
        """
        args: the output list of RAMSES converted by PyMses. 
        output: RAMSES outputs cleanly stored in Zarr format for efficient use.
        """

        root = zarr.open(store_path, mode="w")
        compressor = BloscCodec(
            cname="zstd",
            clevel=3,
            shuffle="bitshuffle",
        )


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
            root.create_array(
                name,
                data=arr,
                chunks=(100_000,),
                compressors=[compressor],
            )

        # --Dust_massdensity-- (even though it is scalar we define a specific chunk because of possible multi-fluid model)
        arr = output["dust_massdensity"]

        if arr.ndim == 1:
            chunks = (50_000,)
        else:
            chunks = (50_000, arr.shape[1])

        root.create_array(
            "dust_massdensity",
            data=arr,
            chunks=chunks,
            compressors=[compressor],
        )

        # --Velocity field--
        if has_velocity:
            arr = output["velocity"]
            root.create_array(
                "velocity",
                data=arr,
                chunks=(50_000, 3),
                compressors=[compressor],
            )

        # --B-fields--
        if has_bfield:
            for component in ("Bx", "By", "Bz"):
                arr = output[component]
                root.create_array(
                    component,
                    data=arr,
                    chunks=(100_000,),
                    compressors=[compressor],
                )


        if has_fluids:
            print('add here the grain sizes and dust number densities')
            arr = output["dust_massdensities"]
            root.create_array(
                "dust_massdensities",
                data=arr,
                chunks=(50_000, arr.shape[1]),
                compressors=[compressor],
            )

        # --Fluid velocity field--
        if has_fluid_v:
            arr = output["fluid_v"]
            root.create_array(
                "fluid_velocity",
                data=arr,
                chunks=(50_000, 3),
                compressors=[compressor],
            )


        # --Dust ratios--
        if has_ratio:
            arr = output["dust_ratio"]
            root.create_array(
                "dust_ratio",
                data=arr,
                chunks=(50_000, arr.shape[1]),
                compressors=[compressor],
            )

        #root.attrs["n_cells"] = N

        root.attrs.update({
            "schema_version": int(1),
            "nb_cells": int(N),
            "nb_species": int(nb_species),
            "l_cm": float(self.l_cm),
            "l_m": float(self.unit_l),
        })

        if ramses_num is not None:
            root.attrs["ramses_output_num"] = int(ramses_num)

        return root


    def ramses(self, ramses_dir, ramses_num, ramses_out, source, sim_param, nb_species):
        importlib.reload(pymses) # Reload pymses to clear internal caches. Another problem with Pymses.

        print(f"\n=== New RAMSES Conversion ===")
        print(f"RAMSES Folder: {ramses_dir}")
        print(f"Output: {ramses_num}")
        print(f"Dust species detected: {nb_species}")
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

        # if nb_species > 0 and nb_species != snap.info["ndust"]:
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
        #     dustratios = np.zeros((nb_species, nr_of_cells))
        #     for i in range(1,nb_species+1):
        #         output["dust_enrich{:d}".format(i)]  = cells["dust_ratio_{:d}".format(i)]  

        #     epsilon_tot = np.zeros(output["dust_enrich1"].shape)
        #     for i in range(1,nb_species+1):
        #         epsilon_tot+=output["dust_enrich{:d}".format(i)]
        #     correction_factor = 1-epsilon_tot
        #     output["gas_massdensity"] = mp*correction_factor*output["density"] # actual gas density 
        
        #     for i in range(1,nb_species+1):
        #         #dustratios[i-1,:] = output["dust_enrich{:d}".format(i)]/correction_factor
        #         output["dust_ratio{:d}".format(i)] = output["dust_enrich{:d}".format(i)]/correction_factor

        #---DENSITY SECTION:---
        #output["density"] = cells["density"]


        if has_ratio:
            #IMPORTANT: Here, ramses dumps the total density = rho_g + rho_d in units of mu micron.(so density is not the actual gas density)
            #for the dust, RAMSES provides dust enrichment (not dust-to-gas mass ratio). So we derive a corretion factor to obtain the dust mass ratios.

            # shape: (cells, species)
            dust_ratio = np.zeros((nr_of_cells, nb_species), dtype=np.float32)
            dust_massdensity = np.zeros((nr_of_cells, nb_species), dtype=np.float32)

            # RAMSES gives enrichment; convert to dust/gas
            epsilon_tot = np.zeros(nr_of_cells, dtype=np.float32)

            enrich = []

            for i in range(1, nb_species + 1):
                e = cells[f"dust_ratio_{i}"]
                enrich.append(e)
                epsilon_tot += e

            correction_factor = 1.0 - epsilon_tot
            output["gas_massdensity"] = mp * correction_factor * cells["density"]

            for i, e in enumerate(enrich):
                dust_ratio[:, i] = e / correction_factor


            for i in range(0, nb_species):
                dust_massdensity[:, i] = dust_ratio[:, i]*output["gas_massdensity"]

            output['dust_massdensities'] = dust_massdensity
            output['dust_massdensity']  = np.sum(dust_massdensity, axis=1)
            
        if has_fluids:
            fluid_density = np.zeros((nr_of_cells, nb_species), dtype=np.float32)

            for i in range(0, nb_species):
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

        #---BFIELD SECTION:---
        if has_bfield:
            # RAMSES uses Constrained Transport: B stored on left and right
            # cell faces. Average the two to get a cell-centred value.
            # Unit conversion: RAMSES code units → Gauss (CGS), as expected by POLARIS.
            if "unit_mag" in snap.info:
                unit_mag = snap.info["unit_mag"].express(C.Gauss)  # T → Gauss
            else:
                # Fallback: derive from standard RAMSES units
                # B_code = sqrt(4π * rho_code * v_code²)  [Gaussian CGS]
                unit_density_cgs = snap.info["unit_density"].express(C.g_cc)
                unit_velocity_cgs = snap.info["unit_velocity"].express(C.cm / C.s)
                unit_mag = np.sqrt(4.0 * np.pi * unit_density_cgs) * unit_velocity_cgs

            B_left  = cells["B_left"]  * unit_mag   # shape (N, 3)
            B_right = cells["B_right"] * unit_mag   # shape (N, 3)
            B_cell  = 0.5 * (B_left + B_right)
            output["Bx"] = B_cell[:, 0]
            output["By"] = B_cell[:, 1]
            output["Bz"] = B_cell[:, 2]

        # #---FLUID VELOCITY SECTION:---
        # if has_fluid_v:
        #     output['fluid_v'] = fluid_cells['fluid_v']

        #---TEMPERATURE SECTION:---
        if has_temp:
            print("\n*************************************************")
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
            scale_l  = pc2cm
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
            print("*************************************************\n")


        root = self.write_amr_to_zarr(output, nb_species,
                                 store_path = ramses_out / f"output{ramses_num}_grid.zarr", 
                                 ramses_num=ramses_num, 
                                 has_temp=has_temp, 
                                 has_velocity=has_velocity, 
                                 has_ratio=has_ratio,
                                 has_fluids=has_fluids,
                                 has_fluid_v=has_fluid_v,
                                 has_bfield=has_bfield)


        return root # will later be stored in the main Grid, then used to write RT files.
        

    def to_radmc(self, ramses_dir, radmc_dir, root, hole_au, f_acc, gridstyle="octtree", coordsystem="cartesian"):
            """
            Convert Zarr grid data to RADMC-3D format.
            
            Args:
                ramses_dir: Path to RAMSES output directory (for sink info)
                radmc_dir: Output directory for RADMC-3D files
                root: Zarr root containing grid data
                hole_au: Radius of hole around sinks (in AU)
                f_acc: Accretion efficiency (unused here, kept for API consistency)
                gridstyle: Grid style for RADMC-3D
                coordsystem: Coordinate system for RADMC-3D
            """
            CLR_LINE = " " * 80 + "\r"
            
            # Read root section
            l_m = root.attrs.get("l_m")
            l_cm = root.attrs.get("l_cm")
            nb_cells = root.attrs.get("nb_cells")
            nb_species = root.attrs.get("nb_species", 1)
            
            level = np.array(root["level"])
            x = np.array(root["x"])
            y = np.array(root["y"])
            z = np.array(root["z"])
            
            # Handle dust density (single species or multi-species)
            if "dust_massdensities" in root:
                dust_massdensity = np.array(root["dust_massdensities"])  # shape: (nb_cells, nb_species)
            elif "dust_massdensity" in root: 
                dust_massdensity = np.array(root["dust_massdensity"])[:, np.newaxis]  # shape: (nb_cells, 1)

            if "gas_massdensity" in root:
                gas_massdensity = np.array(root["gas_massdensity"])[:, np.newaxis]  # shape: (nb_cells, 1)
            
            # Grid bounds (centered at origin)
            x_min = -0.5 * l_m
            y_min = -0.5 * l_m
            z_min = -0.5 * l_m
            
            max_level = int(level.max())
            min_level = int(level.min())
            
            print("\nOctree parameters:")
            print(f"    Level (min, max)     : {min_level}, {max_level}")
            print(f"    Nr. of cells (actual): {nb_cells}")
            print(f"    Nr. of cells (max)   : {8**max_level}")
            print(f"    Length (min, max)    : {l_m/(2**max_level):.3e}, {l_m:.3e} m")
            print(f"    Dust species         : {nb_species}\n")
            
            # Initialize octree with total cell count
            tree = OcTree(x_min, y_min, z_min, l_m)
            tree.nr_of_cells = nb_cells

            # ----------------------------------------------------------
            # Read sinks and dig holes (vectorized)
            # ----------------------------------------------------------
            sinks = ramses.read.sink_info(ramses_dir)

            info_files = sorted(Path(ramses_dir).glob("info_*.txt"))
            for line in open(info_files[0]):
                if "boxlen" in line and "=" in line:
                    boxlen = float(line.split("=")[1])
                    break

            if sinks.num_sinks > 0:
                print(f"Found {sinks.num_sinks} sink(s)")
                x_idx = sinks.columns.index("x")
                y_idx = sinks.columns.index("y")
                z_idx = sinks.columns.index("z")
                sink_positions_m = (sinks.data[:, [x_idx, y_idx, z_idx]] / boxlen - 0.5) * l_m
            else:
                print("No sinks found")

            if sinks.num_sinks > 0 and hole_au > 0:
                hole_radius2 = (hole_au * au2m)**2
                cx = x - 0.5 * l_m
                cy = y - 0.5 * l_m
                cz = z - 0.5 * l_m

                hole_mask = np.zeros(nb_cells, dtype=bool)
                for sp in sink_positions_m:
                    d2 = (cx - sp[0])**2 + (cy - sp[1])**2 + (cz - sp[2])**2
                    hole_mask |= (d2 <= hole_radius2)

                dust_massdensity[hole_mask, :] = 0.0
                print(f"Hole digging: zeroed {hole_mask.sum()} cells "
                    f"around {len(sink_positions_m)} sink(s)")

            # ----------------------------------------------------------
            # Build octree
            # ----------------------------------------------------------
            print("Constructing octree...")
            for i in range(nb_cells):
                c_x = x[i] - 0.5 * l_m
                c_y = y[i] - 0.5 * l_m
                c_z = z[i] - 0.5 * l_m

                cell = CellOct(c_x, c_y, c_z, 0, level[i])
                cell.data = dust_massdensity[i, :].tolist()

                tree.insertInTree(tree.root, cell, 0)

                if i % 10000 == 0:
                    progress = 100.0 * i / nb_cells
                    sys.stdout.write(f"Constructing octree: {progress:.1f}%\r")
                    sys.stdout.flush()

            sys.stdout.write(CLR_LINE)
            print("Constructing octree: done\n")
            
            # Check octree integrity
            print("Checking octree integrity...")
            tree.reset_counter()
            check = tree.checkOcTree(tree.root)
            sys.stdout.write(CLR_LINE)
            
            if not check:
                raise RuntimeError("ERROR: Octree integrity check failed!")
            
            print("Octree structure: OK\n")
            
            # Write octree to RADMC-3D format
            print("Converting to RADMC-3D format...")
            tree.cell_counter = 0
            grid = []
            density = []
            tree._n_species = nb_species
            print('write OcTree call: ')
            tree.writeOcTree_radmc(tree.root, grid, density)
            densityarray = np.array(density)
            sys.stdout.write(CLR_LINE)
            print("Converting to RADMC-3D format: done\n")

            print("Writing the amr_grid.inp file for RADMC-3D...\n")
            radmc3d.write.amr_grid(radmc_dir, 
                                grid, 
                                max_level,
                                nb_cells, 
                                l_cm, 
                                gridstyle=gridstyle, 
                                coordsystem=coordsystem, 
                                x=None, 
                                y=None, 
                                z=None)
            print("Writing amr_grid.inp file: done\n")

            print("Writing the dust_density.inp file for RADMC-3D...\n")
            radmc3d.write.dust_density(radmc_dir, 
                                    densityarray, 
                                    nb_cells,
                                    nb_species, 
                                    gridstyle=gridstyle) 

            print("Writing dust_density.inp file: done\n")


            return grid, densityarray


    def to_polaris(self, ramses_dir, polaris_dir, root, hole_au=0, f_acc=0.1,
                    has_dust_in_sim=True):
            """
            Convert Zarr grid data to a POLARIS binary octree grid file.

            The POLARIS grid stores per-cell data in the following order:
                Bx, By, Bz, Vx, Vy, Vz, gas_mass_density, Tgas, dust_mass_density_1, ..., dust_mass_density_N

            The corresponding POLARIS data IDs written in the header are:
                4 (Bx), 5 (By), 6 (Bz), 7 (Vx), 8 (Vy), 9 (Vz),
                28 (gas density), 3 (gas temperature),
                29 (dust density) x N_dust

            Parameters
            ----------
            ramses_dir : str or Path
                Path to the RAMSES output directory (used to read sink info).
            polaris_dir : str or Path
                Output directory where the POLARIS grid file will be written.
            root : zarr.Group
                Zarr root containing the AMR grid data (from Convert.ramses()).
            hole_au : float
                Radius of the hole around each sink (AU). Density set to 0 inside.
            f_acc : float
                Accretion efficiency factor (not used directly here but kept for
                consistency with to_radmc).
            has_dust_in_sim : bool
                If True, the simulation carries explicit dust fields.
                If False, a single virtual dust species at 1 percent of gas density is assumed.

            Returns
            -------
            output_file : Path
                Path to the written POLARIS binary grid file.
            """
            import struct

            CLR_LINE = " " * 80 + "\r"
            POLARIS_GRID_ID = 20  # octree

            polaris_dir = Path(polaris_dir)
            polaris_dir.mkdir(parents=True, exist_ok=True)

            # ----------------------------------------------------------
            # 0) Read Zarr arrays
            # ----------------------------------------------------------
            l_m = root.attrs.get("l_m")
            l_cm = root.attrs.get("l_cm")
            nb_cells = root.attrs.get("nb_cells")
            nb_species = root.attrs.get("nb_species", 1)

            level = np.array(root["level"])
            x = np.array(root["x"])
            y = np.array(root["y"])
            z = np.array(root["z"])

            # Gas density: Zarr stores g/cc (CGS), POLARIS expects kg/m^3 (SI)
            if "gas_massdensity" in root:
                gas_massdensity = np.array(root["gas_massdensity"]) * 1e3  # g/cc -> kg/m^3
            else:
                raise RuntimeError("gas_massdensity not found in Zarr. "
                                "Ensure the RAMSES grid was created with density fields enabled.")

            # Dust densities: same conversion g/cc -> kg/m^3
            if "dust_massdensities" in root:
                dust_massdensity = np.array(root["dust_massdensities"]) * 1e3
            elif "dust_massdensity" in root:
                dust_massdensity = np.array(root["dust_massdensity"]) * 1e3
                if dust_massdensity.ndim == 1:
                    dust_massdensity = dust_massdensity[:, np.newaxis]
            else:
                raise RuntimeError("No dust density field found in Zarr.")

            n_dust = dust_massdensity.shape[1]

            # B-field: POLARIS expects Gauss (CGS).
            has_bfield = "B_left" in root or "B_right" in root or "Bx" in root
            if has_bfield:
                if "Bx" in root:
                    Bx = np.array(root["Bx"])
                    By = np.array(root["By"])
                    Bz = np.array(root["Bz"])
                else:
                    print("WARNING: B-field arrays not found in expected format. Setting B=0.")
                    has_bfield = False

            if not has_bfield:
                Bx = np.zeros(nb_cells, dtype=np.float32)
                By = np.zeros(nb_cells, dtype=np.float32)
                Bz = np.zeros(nb_cells, dtype=np.float32)

            # Velocity field: Zarr stores m/s (SI), POLARIS expects m/s (SI)
            has_velocity = "velocity" in root
            if has_velocity:
                velocity = np.array(root["velocity"])  # shape (N, 3), already in m/s
                Vx = velocity[:, 0]
                Vy = velocity[:, 1]
                Vz = velocity[:, 2]
            else:
                Vx = np.zeros(nb_cells, dtype=np.float32)
                Vy = np.zeros(nb_cells, dtype=np.float32)
                Vz = np.zeros(nb_cells, dtype=np.float32)

            # Gas temperature (K)
            if "Tgas" in root:
                Tgas = np.array(root["Tgas"])
            elif "Tdust" in root:
                Tgas = np.array(root["Tdust"])
            else:
                print("WARNING: No temperature field found in Zarr. Setting T=10 K everywhere.")
                Tgas = np.full(nb_cells, 10.0, dtype=np.float32)

            # ----------------------------------------------------------
            # 1) Read sinks and dig holes (vectorized)
            # ----------------------------------------------------------
            sinks = ramses.read.sink_info(ramses_dir)

            info_files = sorted(Path(ramses_dir).glob("info_*.txt"))
            for line in open(info_files[0]):
                if "boxlen" in line and "=" in line:
                    boxlen = float(line.split("=")[1])
                    break

            sink_positions_m = None
            if sinks.num_sinks > 0:
                print(f"Found {sinks.num_sinks} sink(s)")
                cols = sinks.columns
                x_col = cols.index("x")
                y_col = cols.index("y")
                z_col = cols.index("z")
                sink_positions_m = (
                    sinks.data[:, [x_col, y_col, z_col]] / boxlen - 0.5
                ) * l_m
            else:
                print("No sinks found")

            if sink_positions_m is not None and hole_au > 0:
                hole_radius2 = (hole_au * au2m) ** 2
                cx = x - 0.5 * l_m
                cy = y - 0.5 * l_m
                cz = z - 0.5 * l_m

                hole_mask = np.zeros(nb_cells, dtype=bool)
                for sp in sink_positions_m:
                    d2 = (cx - sp[0])**2 + (cy - sp[1])**2 + (cz - sp[2])**2
                    hole_mask |= (d2 <= hole_radius2)

                gas_massdensity[hole_mask] = 0.0
                dust_massdensity[hole_mask, :] = 0.0
                print(f"Hole digging: zeroed {hole_mask.sum()} cells "
                    f"around {len(sink_positions_m)} sink(s)")

            # ----------------------------------------------------------
            # 2) Build octree
            # ----------------------------------------------------------
            x_min = -0.5 * l_m
            y_min = -0.5 * l_m
            z_min = -0.5 * l_m

            max_level = int(level.max())
            min_level = int(level.min())

            print("\nPOLARIS Octree parameters:")
            print(f"    Level (min, max)     : {min_level}, {max_level}")
            print(f"    Nr. of cells (actual): {nb_cells}")
            print(f"    Length (min, max)    : {l_m/(2**max_level):.3e}, {l_m:.3e} m")
            print(f"    Dust species         : {n_dust}")
            print(f"    B-field              : {'yes' if has_bfield else 'zeros'}")
            print(f"    Velocity             : {'yes' if has_velocity else 'zeros'}")
            print(f"    Temperature          : {'from Zarr' if 'Tgas' in root else 'default 10 K'}\n")

            tree = OcTree(x_min, y_min, z_min, l_m)
            tree.nr_of_cells = nb_cells

            print("Constructing POLARIS octree...")
            for i in range(nb_cells):
                c_x = x[i] - 0.5 * l_m
                c_y = y[i] - 0.5 * l_m
                c_z = z[i] - 0.5 * l_m

                cell_data = [
                    float(Bx[i]), float(By[i]), float(Bz[i]),
                    float(Vx[i]), float(Vy[i]), float(Vz[i]),
                    float(gas_massdensity[i]),
                    float(Tgas[i]),
                ] + dust_massdensity[i, :].tolist()

                cell = CellOct(c_x, c_y, c_z, 0, int(level[i]))
                cell.data = cell_data

                tree.insertInTree(tree.root, cell, 0)

                if i % 10000 == 0:
                    progress = 100.0 * i / nb_cells
                    sys.stdout.write(f"Constructing POLARIS octree: {progress:.1f}%\r")
                    sys.stdout.flush()

            sys.stdout.write(CLR_LINE)
            print("Constructing POLARIS octree: done\n")

            # Check integrity
            print("Checking octree integrity...")
            tree.reset_counter()
            check = tree.checkOcTree(tree.root)
            sys.stdout.write(CLR_LINE)

            if not check:
                raise RuntimeError("ERROR: POLARIS octree integrity check failed!")
            print("Octree structure: OK\n")

            # ----------------------------------------------------------
            # 3) Write binary POLARIS grid file
            # ----------------------------------------------------------
            ramses_num = int(root.attrs.get("ramses_output_num", 0))
            output_file = polaris_dir / f"ramses_grid_{ramses_num:05d}.dat"

            data_ids = [4, 5, 6, 7, 8, 9, 28, 3] + [29] * n_dust
            data_len = len(data_ids)

            print(f"Writing POLARIS binary grid to: {output_file}")

            with open(output_file, "wb") as f:
                f.write(struct.pack("H", POLARIS_GRID_ID))
                f.write(struct.pack("H", data_len))

                for d_id in data_ids:
                    f.write(struct.pack("H", d_id))

                f.write(struct.pack("d", l_m))

                tree.cell_counter = 0
                tree.writeOcTree(f, tree.root)

            sys.stdout.write(CLR_LINE)
            print("Writing POLARIS grid: done\n")
            print(f"POLARIS octree successfully created: {output_file}\n")

            return output_file


    def to_subboxes(self, ramses_dir, output_dir, root, hole_au, f_acc,
                    which_rad="radmc",
                    box_half_width_au=100.0, isolation_radius_au=100.0,
                    require_luminosity=True, boxlen_pc=None,
                    gridstyle="octtree", coordsystem="cartesian",
                    lmin=1e-3, lmax=1e4, nlam=210, sink_indices=None,
                    min_cells=0):
        """
        Extract one subbox folder per isolated sink, in RADMC-3D or POLARIS format.

        Parameters
        ----------
        ramses_dir : str or Path
            RAMSES output directory (used to read sink info).
        output_dir : str or Path
            Base output directory. Subboxes are written to
            output_dir/subboxes/.
        root : zarr.Group
            Zarr root containing the AMR grid data.
        hole_au : float
            Hole radius around each sink (AU). Density set to 0 inside.
        f_acc : float
            Accretion efficiency factor for computing stellar radii.
        which_rad : str
            Output format: "radmc" for RADMC-3D, "polaris" for POLARIS.
        box_half_width_au : float
            Half-width of each subbox in AU (default 100 = 200 AU boxes).
        isolation_radius_au : float
            Minimum separation for sink filtering in AU.
        require_luminosity : bool
            If True, skip sinks without luminosity data.
        boxlen_pc : float or None
            DEPRECATED — kept for backward compatibility.
            The RAMSES boxlen is now read from the info file automatically.
        gridstyle : str
            "octtree" (default) or "regular". Only used for RADMC-3D.
        coordsystem : str
            "cartesian" (default). Only used for RADMC-3D.
        lmin, lmax : float
            Wavelength range in microns for stars.inp. Only used for RADMC-3D.
        nlam : int
            Number of wavelength points. Only used for RADMC-3D.
        sink_indices : list of int or None
            If provided, skip filtering and process only these sinks.
        min_cells : int
            Minimum number of cells in the subbox for it to be processed.
            Sinks with fewer cells are skipped (default 0 = no minimum).

        Returns
        -------
        results : dict
            {sink_idx: output} for each processed sink.
            For RADMC-3D: output = (grid, densityarray).
            For POLARIS: output = Path to the binary grid file.
        catalog_path : Path
            Path to the CSV catalog of processed sinks.
        """
        from radprocess.pipeline.Subbox import (
            filter_sinks, build_subbox_radmc, build_subbox_polaris,
        )

        if which_rad not in ("radmc", "polaris"):
            raise ValueError(f"which_rad must be 'radmc' or 'polaris', got '{which_rad}'")

        sinks = ramses.read.sink_info(ramses_dir)

        if sinks.num_sinks == 0:
            raise RuntimeError("No sinks found in this RAMSES output.")

        # Read the RAMSES boxlen from the info file.
        # Sink positions in the CSV are in RAMSES code units [0, boxlen].
        # This must match how cell positions were stored in the Zarr
        # (cells.points * boxlen * unit_l).
        ramses_dir = Path(ramses_dir)
        info_files = sorted(ramses_dir.glob("info_*.txt"))
        if info_files:
            with open(info_files[0], "r") as f:
                for line in f:
                    if "boxlen" in line and "=" in line:
                        boxlen = float(line.split("=")[1])
                        break
                else:
                    raise ValueError(f"'boxlen' not found in {info_files[0]}")
        else:
            raise FileNotFoundError(
                f"No info_*.txt found in {ramses_dir}. "
                f"Cannot determine RAMSES boxlen for sink position conversion."
            )

        l_m = root.attrs.get("l_m")

        print(f"RAMSES boxlen: {boxlen:.10f}")
        print(f"Domain size:   {l_m:.6e} m = {l_m/au2m:.1f} AU")
        print(f"Total sinks in output: {sinks.num_sinks}")
        print(f"Output format: {which_rad.upper()}")

        # Filter sinks
        if sink_indices is not None:
            keep = np.array(sink_indices, dtype=int)
            print(f"Using user-provided sink list: {len(keep)} sinks")
        else:
            keep = filter_sinks(
                sinks,
                isolation_radius_au=isolation_radius_au,
                require_luminosity=require_luminosity,
            )
            print(f"After filtering (isolation={isolation_radius_au} AU, "
                  f"require_lum={require_luminosity}): {len(keep)} sinks retained")

        if len(keep) == 0:
            raise RuntimeError("No sinks survived the filtering step.")

        # Prepare output structure
        import csv

        base_dir = Path(output_dir) / "subboxes"
        base_dir.mkdir(parents=True, exist_ok=True)

        cols = sinks.columns
        id_col = cols.index("Id")
        x_col = cols.index("x")
        y_col = cols.index("y")
        z_col = cols.index("z")
        m_col = cols.index("M[Msol]")
        acclum_col = cols.index("acc_lum[Lsol]")
        intlum_col = cols.index("int_lum[Lsol]")

        catalog_path = base_dir / "sink_catalog.csv"
        with open(catalog_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "folder", "sink_id", "row_idx",
                "x_pc", "y_pc", "z_pc",
                "mass_msun", "acc_lum_lsun", "int_lum_lsun",
            ])
            for idx in keep:
                sink_id = int(sinks.data[idx, id_col])
                folder_name = f"sink_{sink_id:04d}"
                writer.writerow([
                    folder_name, sink_id, idx,
                    sinks.data[idx, x_col],
                    sinks.data[idx, y_col],
                    sinks.data[idx, z_col],
                    sinks.data[idx, m_col],
                    sinks.data[idx, acclum_col],
                    sinks.data[idx, intlum_col],
                ])

        print(f"Sink catalog written to: {catalog_path}")

        # Loop over sinks
        results = {}
        skipped = 0

        # Pre-compute centered coordinates and mask components for min_cells check
        l_m = root.attrs.get("l_m")
        x_all = np.array(root["x"]) - 0.5 * l_m
        y_all = np.array(root["y"]) - 0.5 * l_m
        z_all = np.array(root["z"]) - 0.5 * l_m

        x_col = cols.index("x")
        y_col = cols.index("y")
        z_col = cols.index("z")
        hw_m = box_half_width_au * au2m

        for i, idx in enumerate(keep):
            sink_id = int(sinks.data[idx, id_col])

            # Quick cell count check before expensive octree build
            if min_cells > 0:
                sink_pos_m = (sinks.data[idx, [x_col, y_col, z_col]] / boxlen - 0.5) * l_m
                n_cells = int(np.sum(
                    (np.abs(x_all - sink_pos_m[0]) <= hw_m) &
                    (np.abs(y_all - sink_pos_m[1]) <= hw_m) &
                    (np.abs(z_all - sink_pos_m[2]) <= hw_m)
                ))
                if n_cells < min_cells:
                    print(f"\n  Skipping sink {sink_id}: only {n_cells} cells "
                          f"(min_cells={min_cells})")
                    skipped += 1
                    continue

            print(f"\n{'='*60}")
            print(f"  Processing sink {sink_id}  ({i+1-skipped}/{len(keep)-skipped})")
            print(f"{'='*60}")

            sink_dir = base_dir / f"sink_{sink_id:04d}"

            if which_rad == "radmc":
                result = build_subbox_radmc(
                    root=root,
                    ramses_dir=ramses_dir,
                    output_dir=sink_dir,
                    sink_idx=idx,
                    sinks=sinks,
                    box_half_width_au=box_half_width_au,
                    hole_au=hole_au,
                    f_acc=f_acc,
                    boxlen=boxlen,
                    gridstyle=gridstyle,
                    coordsystem=coordsystem,
                    lmin=lmin,
                    lmax=lmax,
                    nlam=nlam,
                )
            else:
                result = build_subbox_polaris(
                    root=root,
                    ramses_dir=ramses_dir,
                    output_dir=sink_dir,
                    sink_idx=idx,
                    sinks=sinks,
                    box_half_width_au=box_half_width_au,
                    hole_au=hole_au,
                    f_acc=f_acc,
                    boxlen=boxlen,
                )

            results[sink_id] = result

        print(f"\n{'='*60}")
        print(f"  All done: {len(results)} subboxes written to {base_dir}")
        if skipped > 0:
            print(f"  Skipped: {skipped} sinks (fewer than {min_cells} cells)")
        print(f"{'='*60}\n")

        return results, catalog_path


    def polaris_photon(self):
        self.rat1 = 1
