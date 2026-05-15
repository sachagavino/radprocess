import os
import json

from dataclasses import fields
from pathlib import Path

import zarr

from radprocess.pipeline.Grid import Grid
from radprocess.pipeline.Convert import Convert
from radprocess import radmc3d
from radprocess import ramses
from radprocess.utils.config import ConfigParams

class Pipeline:

    def __init__(self):
        self.convert = Convert() 
        self.grid = Grid()
        self.configparams = ConfigParams()

    @property
    def ramses_converted_dir(self):
        d = Path(self.configparams.dir.pipeline_output) / "ramses"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def radmc_outputs_dir(self):
        d = Path(self.configparams.dir.pipeline_output) / "radmc3d"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def polaris_outputs_dir(self):
        d = Path(self.configparams.dir.pipeline_output) / "polaris"
        d.mkdir(parents=True, exist_ok=True)
        return d

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
        ramses_dir = self.configparams.dir.ramses_output
        try:
            nvar, variables, nb_dust = ramses.read.hydro_file_descriptor(ramses_dir)

            text = f"nvar = {nvar}\n\nVariables:\n"
            for idx, name in variables.items():
                text += f"  #{idx}: {name}\n"

            # text += f"\nDust ratios detected: {nb_dust}\n"
            text += f"\nDust components detected: {nb_dust}\n"
            self.configparams.nb_dust = nb_dust
            return text

        except Exception as e:
            return f"Error reading hydro_file_descriptor.txt:\n{e}"

    def read_other_file_descriptor(self):
        """
        Reads other *_file_descriptor.txt (mf_, mp_, etc.) from RAMSES directory.
        If found, overrides nb_dust using number of fluids.
        """
        ramses_dir = self.configparams.dir.ramses_output

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

        Returns
        -------
        SinkInfo
            A dataclass with columns, rows, num_sinks, and data.
            Renders as a collapsible HTML table in Jupyter notebooks.
        """

        ramses_dir = self.configparams.dir.ramses_output
        try:
            return ramses.read.sink_info(ramses_dir)
        except Exception as e:
            print(f"Error reading sink info: {e}")
            return None
        
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
    #         self.configparams.dir.ramses_output
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
        ramses_dir = self.configparams.dir.ramses_output
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
        self.configparams.dir.pipeline_output = str(base)

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
        ramses_dir = self.configparams.dir.ramses_output
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
        ramses_dir = self.configparams.dir.ramses_output
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

    def convert_to_polaris(self, hole_au=None):
        """
        Convert the RAMSES Zarr grid to a POLARIS binary octree grid.

        Parameters
        ----------
        hole_au : float or None
            Hole radius around sinks (AU). If None, uses configparams.sim.size_hole_au.

        Returns
        -------
        output_file : Path
            Path to the written POLARIS binary grid file.
        """
        ramses_dir = self.configparams.dir.ramses_output
        polaris_dir = self.polaris_outputs_dir
        if hole_au is None:
            hole_au = self.configparams.sim.size_hole_au
        f_acc = self.configparams.sim.facc

        root = self.get_amr_root()

        num = int(root.attrs["ramses_output_num"])
        print(f"Using RAMSES output #{num} from Zarr")

        output_file = self.convert.to_polaris(
            ramses_dir,
            polaris_dir,
            root,
            hole_au=hole_au,
            f_acc=f_acc,
        )

        return output_file

    def convert_subboxes(self, which_rad="radmc",
                         box_half_width_au=100.0, isolation_radius_au=100.0,
                         hole_au=None, boxlen_pc=None, require_luminosity=True,
                         sink_indices=None, gridstyle="octtree", coordsystem="cartesian"):
        """
        Extract one subbox folder per isolated sink, in RADMC-3D or POLARIS format.

        Parameters
        ----------
        which_rad : str
            Output format: "radmc" for RADMC-3D, "polaris" for POLARIS.
        box_half_width_au : float
            Half-width of each subbox in AU (default 100 = 200 AU boxes).
        isolation_radius_au : float
            Minimum separation for sink filtering in AU.
        hole_au : float or None
            Hole radius around each sink (AU). If None, uses configparams.sim.size_hole_au.
        boxlen_pc : float or None
            RAMSES box size in pc. If None, derived from Zarr.
        require_luminosity : bool
            If True, skip sinks without luminosity data.
        sink_indices : list of int or None
            If provided, skip filtering and process only these sinks.
        gridstyle : str
            "octtree" (default) or "regular". Only used for RADMC-3D.
        coordsystem : str
            "cartesian" (default). Only used for RADMC-3D.

        Returns
        -------
        results : dict
            {sink_idx: output} for each processed sink.
        catalog_path : Path
            Path to the CSV catalog of processed sinks.
        """
        ramses_dir = self.configparams.dir.ramses_output
        if hole_au is None:
            hole_au = self.configparams.sim.size_hole_au
        f_acc = self.configparams.sim.facc

        # Choose output directory based on format
        if which_rad == "radmc":
            output_dir = self.radmc_outputs_dir
        elif which_rad == "polaris":
            output_dir = self.polaris_outputs_dir
        else:
            raise ValueError(f"which_rad must be 'radmc' or 'polaris', got '{which_rad}'")

        root = self.get_amr_root()

        return self.convert.to_subboxes(
            ramses_dir=ramses_dir,
            output_dir=output_dir,
            root=root,
            hole_au=hole_au,
            f_acc=f_acc,
            which_rad=which_rad,
            box_half_width_au=box_half_width_au,
            isolation_radius_au=isolation_radius_au,
            require_luminosity=require_luminosity,
            boxlen_pc=boxlen_pc,
            gridstyle=gridstyle,
            coordsystem=coordsystem,
            sink_indices=sink_indices,
        )

    def run_polaris_opacity(
        self,
        dust_components,
        dust_size_min=None,
        dust_size_max=None,
        dust_size_powerlaw=None,
        mean_molecular_weight=None,
        mass_fraction=None,
        nr_threads=None,
        grid_path=None,
        n_dust_override=None,
        polaris_binary=None,
        cleanup=True,
    ):
        """
        Run POLARIS with 1 photon package to generate dust opacity tables.

        This is Step 4 of the pipeline. It requires that a POLARIS grid
        file already exists (from convert_to_polaris, Step 2).

        Parameters default to values in configparams.polaris when not
        explicitly provided. dust_components must always be provided
        explicitly since it contains file paths specific to the user's setup.

        Parameters
        ----------
        dust_components : list of dict
            Dust material definitions. Each dict has:
                path (str): path to the .nk or .cs file
                weight (float): mass fraction weight
            Example:
                [{"path": "/path/silicate.cs", "weight": 0.625},
                 {"path": "/path/carbon.cs",   "weight": 0.375}]
        dust_size_min : float or None
            Minimum grain radius in metres.
        dust_size_max : float or None
            Maximum grain radius in metres.
        dust_size_powerlaw : float or None
            Power-law exponent (default -3.5 for MRN).
        mean_molecular_weight : float or None
            Gas mean molecular weight.
        mass_fraction : float or None
            Dust-to-gas mass fraction.
        nr_threads : int or None
            Number of OpenMP threads.
        grid_path : str or Path or None
            Path to the POLARIS grid file. If None, auto-detected.
        n_dust_override : int or None
            Override the number of dust species.
        polaris_binary : str or None
            Name or path of the POLARIS executable.
        cleanup : bool
            If True, remove previous POLARIS run outputs before starting.

        Returns
        -------
        data_dir : Path
            Path to the POLARIS data/ directory containing dust_mixture_*.dat.
        """
        from radprocess.polaris.opacity import run_opacity

        ramses_dir = self.configparams.dir.ramses_output
        polaris_dir = self.polaris_outputs_dir
        f_acc = self.configparams.sim.facc
        pc = self.configparams.polaris

        # Use dataclass defaults for any parameter not explicitly set
        if dust_size_min is None:
            dust_size_min = pc.dust_size_min
        if dust_size_max is None:
            dust_size_max = pc.dust_size_max
        if dust_size_powerlaw is None:
            dust_size_powerlaw = pc.dust_size_powerlaw
        if mean_molecular_weight is None:
            mean_molecular_weight = pc.mean_molecular_weight
        if mass_fraction is None:
            mass_fraction = pc.mass_fraction
        if nr_threads is None:
            nr_threads = pc.nr_threads
        if polaris_binary is None:
            polaris_binary = pc.polaris_binary

        # Auto-detect grid file if not provided
        if grid_path is None:
            candidates = list(polaris_dir.glob("ramses_grid_*.dat"))
            if not candidates:
                raise FileNotFoundError(
                    f"No POLARIS grid file found in {polaris_dir}. "
                    "Run convert_to_polaris() first (Step 2)."
                )
            grid_path = candidates[0]
            print(f"Auto-detected POLARIS grid: {grid_path}")

        data_dir = run_opacity(
            ramses_dir=ramses_dir,
            polaris_dir=polaris_dir,
            grid_path=grid_path,
            dust_components=dust_components,
            dust_size_min=dust_size_min,
            dust_size_max=dust_size_max,
            dust_size_powerlaw=dust_size_powerlaw,
            mean_molecular_weight=mean_molecular_weight,
            mass_fraction=mass_fraction,
            nr_threads=nr_threads,
            f_acc=f_acc,
            n_dust_override=n_dust_override,
            polaris_binary=polaris_binary,
            cleanup=cleanup,
        )

        return data_dir

    def prepare_radmc3d_inputs(
        self,
        polaris_data_dir=None,
        n_dust=None,
        wave_min=None,
        wave_max=None,
        n_wavelengths=None,
        nphot=None,
        setthreads=None,
        scattering_mode=None,
        polaris_skiprows=29,
    ):
        """
        Convert POLARIS opacities to RADMC-3D format and write all
        remaining RADMC-3D input files (Step 5).

        Requires that Steps 3 (RADMC-3D grid) and 4 (POLARIS opacity run)
        have already been completed.

        Parameters default to values in configparams.radmc3d when not
        explicitly provided.

        Parameters
        ----------
        polaris_data_dir : str or Path or None
            Path to the POLARIS data/ directory. If None, auto-detected
            from {polaris_outputs_dir}/data/.
        n_dust : int or None
            Number of dust species. If None, auto-detected from POLARIS output.
        wave_min, wave_max : float or None
            Wavelength range in microns.
        n_wavelengths : int or None
            Number of wavelength points.
        nphot : int or None
            Number of photon packages for mctherm.
        setthreads : int or None
            Number of OpenMP threads.
        scattering_mode : int or None
            Scattering mode (1 = isotropic).
        polaris_skiprows : int
            Header rows to skip in POLARIS opacity files.

        Returns
        -------
        radmc_dir : Path
            The RADMC-3D directory containing all input files.
        """
        from radprocess.radmc3d.prepare import prepare_radmc3d_inputs

        ramses_dir = self.configparams.dir.ramses_output
        radmc_dir = self.radmc_outputs_dir
        f_acc = self.configparams.sim.facc
        rc = self.configparams.radmc3d

        # Use dataclass defaults for any parameter not explicitly set
        if wave_min is None:
            wave_min = rc.wave_min
        if wave_max is None:
            wave_max = rc.wave_max
        if n_wavelengths is None:
            n_wavelengths = rc.n_wavelengths
        if nphot is None:
            nphot = rc.nphot
        if setthreads is None:
            setthreads = rc.setthreads
        if scattering_mode is None:
            scattering_mode = rc.scattering_mode

        if polaris_data_dir is None:
            polaris_data_dir = self.polaris_outputs_dir / "data"
            if not polaris_data_dir.exists():
                raise FileNotFoundError(
                    f"POLARIS data directory not found: {polaris_data_dir}. "
                    "Run run_polaris_opacity() first (Step 4)."
                )

        return prepare_radmc3d_inputs(
            ramses_dir=ramses_dir,
            radmc_dir=radmc_dir,
            polaris_data_dir=polaris_data_dir,
            n_dust=n_dust,
            f_acc=f_acc,
            wave_min=wave_min,
            wave_max=wave_max,
            n_wavelengths=n_wavelengths,
            nphot=nphot,
            setthreads=setthreads,
            scattering_mode=scattering_mode,
            polaris_skiprows=polaris_skiprows,
        )

    def run_radmc3d_mctherm(self, radmc3d_binary="radmc3d"):
        """
        Execute ``radmc3d mctherm`` to compute the dust temperature (Step 6).

        All input files must already exist in the RADMC-3D output directory
        (from Steps 3 and 5).

        Parameters
        ----------
        radmc3d_binary : str
            Name or path of the RADMC-3D executable.

        Returns
        -------
        temp_file : Path
            Path to the dust_temperature.bdat output file.
        """
        from radprocess.radmc3d.run import mctherm

        radmc_dir = self.radmc_outputs_dir

        # Infer output number for the log filename
        root = self.get_amr_root()
        output_num = int(root.attrs.get("ramses_output_num", 0))
        log_path = radmc_dir / f"radmc3d_mctherm_{output_num:05d}.log"

        temp_file = mctherm(
            radmc_dir=radmc_dir,
            log_path=log_path,
            radmc3d_binary=radmc3d_binary,
        )

        return temp_file

    def merge_temperature(self, n_dust=None):
        """
        Merge RADMC-3D dust temperatures into the POLARIS grid (Step 7).

        Reads the POLARIS grid_temp.dat (from Step 4) and the RADMC-3D
        dust_temperature.bdat (from Step 6), replaces the POLARIS dust
        temperatures with the RADMC-3D values, and writes the result as
        grid_temp.radmc3d.dat in the POLARIS output directory.

        Parameters
        ----------
        n_dust : int or None
            Number of dust species. If None, auto-detected from the
            RAMSES info file.

        Returns
        -------
        merged_grid : Path
            Path to the merged grid file (grid_temp.radmc3d.dat).
        """
        from radprocess.polaris.merge import merge_radmc3d_temperature

        ramses_dir = self.configparams.dir.ramses_output
        polaris_dir = self.polaris_outputs_dir
        radmc_dir = self.radmc_outputs_dir

        return merge_radmc3d_temperature(
            polaris_dir=polaris_dir,
            radmc_dir=radmc_dir,
            n_dust=n_dust,
            ramses_dir=ramses_dir,
        )

    def render_images(
        self,
        dust_components,
        npix,
        distance_pc,
        wavelengths_mm,
        views=None,
        midplane_zoom=1,
        fov_m=None,
        label="whole",
        grid_path=None,
        n_dust=None,
        nr_threads=None,
        dust_size_min=None,
        dust_size_max=None,
        dust_size_powerlaw=None,
        mean_molecular_weight=None,
        mass_fraction=None,
        polaris_binary=None,
        cleanup_views=True,
        polaris_cmd="CMD_DUST_EMISSION",
        alignment="ALIG_PA",
        peel_off=True,
        acceptance_angle=None,
        nr_photons_scat=None,
        source_star_scat=None,
    ):
        """
        Run POLARIS dust emission/scattering imaging (Step 8).

        Renders images at the specified wavelengths and viewing angles
        using the merged POLARIS grid (from Step 7).

        Parameters
        ----------
        dust_components : list of dict
            Dust material definitions (same as for the opacity run).
        npix : int
            Image resolution (npix x npix pixels).
        distance_pc : float
            Source distance in parsecs.
        wavelengths_mm : list of float
            Wavelengths to image in millimetres.
        views : list of str or None
            Viewing angles to render (e.g. ["xy", "xz", "yz"]).
            If None, renders all three standard views.
        midplane_zoom : int or float
            Midplane zoom factor (1 = full box, 10 = zoomed inner).
        fov_m : float or None
            Field of view in metres (half-width). If None, POLARIS uses
            the full grid extent. Set this for zoomed/inner imaging.
        label : str
            Output subdirectory label ("whole" or "inner").
        grid_path : str or Path or None
            Path to the merged grid. If None, auto-detected.
        n_dust : int or None
            Number of dust species. If None, auto-detected from RAMSES info.
        nr_threads : int or None
            OpenMP threads (defaults to configparams.polaris.nr_threads).
        dust_size_min, dust_size_max : float or None
            Grain size range in metres.
        dust_size_powerlaw : float or None
            Size distribution exponent.
        mean_molecular_weight : float or None
            Gas mu.
        mass_fraction : float or None
            Dust-to-gas mass fraction.
        polaris_binary : str or None
            POLARIS executable name or path.
        cleanup_views : bool
            Remove previous image outputs before rendering.
        polaris_cmd : str
            POLARIS command: "CMD_DUST_EMISSION", "CMD_DUST_SCATTERING",
            or both separated by a space.
        alignment : str
            Grain alignment: "ALIG_PA", "ALIG_IDG", "ALIG_RAT",
            "ALIG_INTERNAL", or "" for no alignment.
        peel_off : bool
            Use peel-off technique for scattering.
        acceptance_angle : float or None
            Acceptance angle for scattered light (degrees).
        nr_photons_scat : int or None
            Number of photon packages for scattering Monte Carlo.
        source_star_scat : list of dict or None
            Stellar sources for scattering (pos_m, radius_rsun, temperature_k).

        Returns
        -------
        image_dirs : dict
            {view_name: Path} for each rendered view.
        """
        from radprocess.polaris.imaging import render_images
        from radprocess.polaris.opacity import _find_info_file, _read_ndust

        polaris_dir = self.polaris_outputs_dir
        image_output_dir = Path(self.configparams.dir.pipeline_output) / "images"
        pc = self.configparams.polaris

        # Defaults from config dataclass
        if dust_size_min is None:
            dust_size_min = pc.dust_size_min
        if dust_size_max is None:
            dust_size_max = pc.dust_size_max
        if dust_size_powerlaw is None:
            dust_size_powerlaw = pc.dust_size_powerlaw
        if mean_molecular_weight is None:
            mean_molecular_weight = pc.mean_molecular_weight
        if mass_fraction is None:
            mass_fraction = pc.mass_fraction
        if nr_threads is None:
            nr_threads = pc.nr_threads
        if polaris_binary is None:
            polaris_binary = pc.polaris_binary

        # Auto-detect grid
        if grid_path is None:
            grid_path = polaris_dir / "grid_temp.radmc3d.dat"
            if not grid_path.exists():
                raise FileNotFoundError(
                    f"Merged grid not found: {grid_path}. "
                    "Run merge_temperature() first (Step 7)."
                )

        # Auto-detect n_dust
        if n_dust is None:
            ramses_dir = self.configparams.dir.ramses_output
            info_path = _find_info_file(ramses_dir)
            n_dust = _read_ndust(info_path)
            if n_dust == 0:
                n_dust = 1

        # Get output number for naming
        root = self.get_amr_root()
        output_num = int(root.attrs.get("ramses_output_num", 0))

        return render_images(
            polaris_dir=polaris_dir,
            image_output_dir=image_output_dir,
            grid_path=grid_path,
            dust_components=dust_components,
            n_dust=n_dust,
            dust_size_min=dust_size_min,
            dust_size_max=dust_size_max,
            dust_size_powerlaw=dust_size_powerlaw,
            mean_molecular_weight=mean_molecular_weight,
            mass_fraction=mass_fraction,
            npix=npix,
            distance_pc=distance_pc,
            wavelengths_mm=wavelengths_mm,
            views=views,
            nr_threads=nr_threads,
            midplane_zoom=midplane_zoom,
            fov_m=fov_m,
            output_num=output_num,
            polaris_binary=polaris_binary,
            label=label,
            cleanup_views=cleanup_views,
            polaris_cmd=polaris_cmd,
            alignment=alignment,
            peel_off=peel_off,
            acceptance_angle=acceptance_angle,
            nr_photons_scat=nr_photons_scat,
            source_star_scat=source_star_scat,
        )