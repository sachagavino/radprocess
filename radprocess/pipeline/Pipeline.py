import os
import json

from dataclasses import fields
from pathlib import Path

import zarr

from radprocess.pipeline.Grid import Grid
from radprocess.pipeline.Convert import Convert
from radprocess.pipeline.Subbox import convert_subboxes_to_radmc
from radprocess import radmc3d
from radprocess import ramses
from radprocess.utils.config import ConfigParams

from radprocess.constants.constants import pc2m

class Pipeline:

    def __init__(self):
        self.convert = Convert() 
        self.grid = Grid()
        self.configparams = ConfigParams()

    @property
    def ramses_converted_dir(self):
        return Path(self.configparams.pipoutput.main_output_dir) / "ramses"

    @property
    def radmc_outputs_dir(self):
        return Path(self.configparams.pipoutput.main_output_dir) / "radmc3d"

    @property
    def polaris_outputs_dir(self):
        return Path(self.configparams.pipoutput.main_output_dir) / "polaris"

    def get_amr_root(self):
        """
        Locate and return the RAMSES AMR Zarr root.

        Order:
          1) in-memory Grid
          2) on disk in ramses_converted_dir

        Also validates presence of required metadata.

        Returns
        -------
        root : zarr.Group

        Raises
        ------
        RuntimeError if not found or invalid.
        """

        root = None

        # --------------------------------------------------
        # 1) In-memory
        # --------------------------------------------------
        if self.grid.amr_grid:
            root = self.grid.amr_grid[0]

        # --------------------------------------------------
        # 2) Disk
        # --------------------------------------------------
        if root is None:

            zarr_dir = Path(self.ramses_converted_dir)

            if not zarr_dir.exists():
                raise RuntimeError(
                    "No RAMSES Zarr directory found. "
                    "Run 'Create RAMSES grid' first."
                )

            candidates = list(zarr_dir.glob("*.zarr"))

            if not candidates:
                raise RuntimeError(
                    "No Zarr files found in ramses_converted_dir."
                )

            # For now assume single Zarr
            root = zarr.open(candidates[0], mode="r")

            # Cache for later steps
            self.grid.amr_grid = [root]

        # --------------------------------------------------
        # 3) Validate
        # --------------------------------------------------
        if "ramses_output_num" not in root.attrs:
            raise RuntimeError(
                "Zarr root missing 'ramses_output_num' attribute."
            )

        if "nb_cells" not in root.attrs:
            raise RuntimeError(
                "Zarr root missing 'nb_cells' attribute."
            )

        return root

    def get_amr_fields(self):
        """
        Return available scalar/vector fields from AMR Zarr root.
        """

        root = self.get_amr_root()   # <-- IMPORTANT (loads from disk if needed)

        # root is now guaranteed to exist or raise

        fields = []

        for key in root.array_keys():
            # skip coordinates if you want
            if key in ("x", "y", "z", "dx", "level"):
                continue
            fields.append(key)

        return sorted(fields)



    def get_enabled_amr_fields(self):

        if getattr(self.configparams.amrsource, "temp", False):
            if not getattr(self.configparams.amrsource, "p", False):
                print("NOTE: 'temp' requires 'pressure'")
                self.configparams.amrsource.p = True

        enabled_vars = []
        dust_count = self.configparams.nb_dust

        for f in fields(self.configparams.amrsource):
            if not getattr(self.configparams.amrsource, f.name):
                continue

            meta = f.metadata or {}
            ramses_name = meta.get("ramses_name", f.name)
            field_type = meta.get("type", "scalar")

            # ---------- dust ----------
            if ramses_name == "dust_ratio":
                for i in range(1, dust_count + 1):
                    enabled_vars.append(f"dust_ratio_{i}")
                continue

            # ---------- fluid density ----------
            if ramses_name == "fluid_density":
                for i in range(1, dust_count + 1):
                    enabled_vars.append(f"fluid_density_{i}")
                continue

            # ---------- fluid velocity (vector) ----------
            if ramses_name == "fluid_v":
                for i in range(1, dust_count + 1):
                    enabled_vars.append(f"fluid_v_{i}")
                continue

            # ---------- normal vector ----------
            if field_type == "vector":
                enabled_vars.append(ramses_name)
                continue

            # ---------- scalar ----------
            enabled_vars.append(ramses_name)

        return enabled_vars


    def read_hydro_descriptor(self):
        """
        Reads hydro_file_descriptor.txt from the current RAMSES directory.
        Returns a formatted string.
        """
        ramses_dir = self.configparams.ramsesoutput.ramses_output_dir
        try:
            nvar, variables, nb_dust = ramses.read.hydro_file_descriptor(ramses_dir)

            text = f"nvar = {nvar}\n\nVariables:\n"
            for idx, name in variables.items():
                text += f"  #{idx}: {name}\n"

            text += f"\nDust ratios detected: {nb_dust}\n"
            self.configparams.nb_dust = nb_dust
            return text

        except Exception as e:
            return f"Error reading hydro_file_descriptor.txt:\n{e}"

    def read_other_file_descriptor(self):
        """
        Reads other *_file_descriptor.txt (mf_, mp_, etc.) from RAMSES directory.
        If found, overrides nb_dust using number of fluids.
        """
        ramses_dir = self.configparams.ramsesoutput.ramses_output_dir

        try:
            nvar, variables, nb_fluids = ramses.read.other_file_descriptor(ramses_dir)

            if nvar is None:
                return "No additional *_file_descriptor.txt found."

            text = f"nvar = {nvar}\n\nVariables:\n"
            for idx, name in variables.items():
                text += f"  #{idx}: {name}\n"

            text += f"\nFluids detected: {nb_fluids}\n"

            # IMPORTANT: override hydro result
            self.configparams.nb_dust = nb_fluids

            return text

        except Exception as e:
            return f"Error reading other file_descriptor.txt:\n{e}"

        

    def read_sink_info(self):
        """
        Reads the sink_0000.info file from the RAMSES output directory
        using ramses.read.sink_info() function.
        Returns a string suitable for HTML display.
        """

        ramses_dir = self.configparams.ramsesoutput.ramses_output_dir
        try:
            sink_data = ramses.read.sink_info(ramses_dir)
        except Exception as e:
            return f"Error reading sink info: {e}"

        # Convert to aligned string
        columns = sink_data.columns
        rows = sink_data.rows

        lines = ["  ".join(columns)]
        for row in rows:
            lines.append("  ".join(str(row[col]) for col in columns))

        return "\n".join(lines)
        
    def read_pymsesrc(self):
        """Reads the ~/.pymses/pymsesrc JSON file and returns a dictionary."""
        pymses_directory = os.path.expanduser("~/.pymses")
        filename = os.path.join(pymses_directory, "pymsesrc")

        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} not found. Please create the pymsesrc file.")

        try:
            with open(filename, "r") as f:
                config_dict = json.load(f)
        except json.JSONDecodeError as e:
            raise RuntimeError(f"Error decoding JSON in {filename}:\n{e}")
        except Exception as e:
            raise RuntimeError(f"Error while reading {filename}:\n{e}")

        return config_dict

    # def set_pymsesrc(self):
    #     """
    #     Write ~/.pymses/pymsesrc based on current configuration,
    #     then return the Pymsesrc object so Jupyter displays it.
    #     """
    #     self.ndust = ramses.read.hydro_file_descriptor(
    #         self.configparams.ramsesoutput.ramses_output_dir
    #     )[2]
    #     print(f'there is {ndust} dust species in the RAMSES simulation.')  
    #     # Write the file
    #     ramses.write.pymsesrc(self,
    #         ndust=self.ndust,
    #         rho=True,
    #         dustratios=self.configparams.pymsesrc.dustratios,
    #         vel=self.configparams.pymsesrc.vel,
    #         bl=self.configparams.pymsesrc.bl,
    #         br=self.configparams.pymsesrc.br,
    #         p=self.configparams.pymsesrc.p,
    #         xi=self.configparams.pymsesrc.xi,
    #         phi=self.configparams.pymsesrc.phi,
    #         g=self.configparams.pymsesrc.g,
    #     )

    #     # RETURN Pymsesrc config block → Jupyter displays it
    #     return self.configparams.pymsesrc

    def set_pymsesrc(self):
        """
        Write ~/.pymses/pymsesrc based on current configuration.
        """
        ramses_dir = self.configparams.ramsesoutput.ramses_output_dir
        #ndust = ramses.read.hydro_file_descriptor(ramses_dir)[2]
        #print(f"There are {ndust} dust species in this RAMSES simulation.")
        # Write the file
        ramses.write.pymsesrc(ramses_dir)
    
    def set_working_dir(self, workdir):
        """
        Set main pipeline output directory and create subfolders:
          - ramses
          - radmc3d
          - polaris
        """

        base = Path(workdir).expanduser().resolve()

        # Create main dir
        base.mkdir(parents=True, exist_ok=True)

        # Subdirectories
        for sub in ["ramses", "radmc3d", "polaris"]:
            (base / sub).mkdir(exist_ok=True)

        # Store ONLY the main dir in config
        self.configparams.pipoutput.main_output_dir = str(base)

        return (
            f"✔ Working directory set to:\n{base}\n\n"
            f"Created/verified:\n"
            f" - {base / 'ramses'}\n"
            f" - {base / 'radmc3d'}\n"
            f" - {base / 'polaris'}"
        )

  
    def thermal_radmc3d(self,
                        run=True, 
                        nphot=1e4, 
                        write_opac=True,
                        write_control=True, 
                        write_star=True, 
                        write_wave=True, 
                        write_mcmono=True, 
                        write_ext=True, 
                        **keywords):
        """ 
        Notes:
        run MC dust radiative transfer, open the resulting dust temperature as an array 
        and computes the surface-area weigthed temperature. If run == False, user assumes
        the RADMC3D output files already exist.
        -----
	    """	
        self.write_radmc3d(nphot_therm=nphot, 
                           write_opac=write_opac, 
                           write_control=write_control, 
                           write_star=write_star, 
                           write_wave=write_wave, 
                           write_mcmono=write_mcmono, 
                           write_ext=write_ext, 
                           **keywords)

        if write_control == False or write_opac == False or write_star == False or write_wave == False:
            print('WARNING: most RADMC3D input files will not be created. Will continue..\
                   but errors can be raised if one or more required input files are missing.\n')

        if run == True:
            self.run_thermal_radmc3d(nphot=nphot, **keywords)


    def write_radmc3d(self, 
                      write_dens, 
                      write_grid, 
                      write_opac, 
                      write_control, 
                      write_star, 
                      write_wave, 
                      write_mcmono, 
                      write_ext, 
                      **keywords):
        print('\nWRITING RADMC3D INPUT FILES:')
        print('----------------------------\n')
        #os.system("rm thermal/*.inp")
        if not os.path.exists('thermal'):
            os.makedirs('thermal')

        if write_control==True:
            radmc3d.write.control(**keywords)

    def load_ramses(self):
        """
        read RAMSES files, convert into AMR and store
        the grid inside the amr_grid from Grid.
        """
        ramses_dir = self.configparams.ramsesoutput.ramses_output_dir
        #radmc_dir = self.configparams.pipoutput.radmc_output_dir
        ramses_out = self.ramses_converted_dir
        ramses_out.mkdir(parents=True, exist_ok=True)

        enabled_source = self.get_enabled_amr_fields()
        nb_sizes = self.configparams.nb_dust
        sim_param = self.configparams.sim

        clean = ramses_dir.strip().rstrip("/")   # remove whitespace + trailing slash
        folder = clean.rsplit("/", 1)[0] + "/"
        num = int(clean.rsplit("_", 1)[1])

        amr_grid = self.convert.ramses(
            folder,
            num,
            ramses_out,
            enabled_source,
            sim_param,
            nb_sizes,
        )

        self.grid.add_amr_grid(amr_grid)

    def convert_to_radmc(self, gridstyle="octtree", coordsystem="cartesian"):
        ramses_dir = self.configparams.ramsesoutput.ramses_output_dir
        radmc_dir = self.radmc_outputs_dir
        hole_au = self.configparams.sim.size_hole_au
        f_acc = self.configparams.sim.facc

        # --------------------------------------------------
        # 1) Extract RAMSES output number from Zarr
        # --------------------------------------------------
        root = self.get_amr_root()

        num = int(root.attrs["ramses_output_num"])
        print(f"Using RAMSES output #{num} from Zarr")

        # --------------------------------------------------
        # 2) Run and store conversion
        # --------------------------------------------------
        radmc_grid, radmc_dens = self.convert.to_radmc(ramses_dir, 
                                                       radmc_dir, 
                                                       root, 
                                                       hole_au,
                                                       f_acc,  
                                                       gridstyle=gridstyle, 
                                                       coordsystem=coordsystem)
        self.grid.add_radmc_grid(radmc_grid)
        # nb_sizes = self.configparams.nb_dust
        # sim_param = self.configparams.sim

    def convert_subboxes(self, box_half_width_au=100.0, isolation_radius_au=100.0,
                        hole_au=None, boxlen_pc=None, require_luminosity=True):
        if hole_au is None:
            hole_au = self.configparams.sim.size_hole_au
        if boxlen_pc is None:
            root = self.get_amr_root()
            boxlen_pc = root.attrs.get("l_m") / pc2m

        return convert_subboxes_to_radmc(
            self,
            box_half_width_au=box_half_width_au,
            isolation_radius_au=isolation_radius_au,
            hole_au=hole_au,
            boxlen_pc=boxlen_pc,
            require_luminosity=require_luminosity,
        )




