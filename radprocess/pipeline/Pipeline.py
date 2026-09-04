import os
import json

from dataclasses import fields
from pathlib import Path

import zarr
import numpy as np

from radprocess.pipeline.Grid import Grid
from radprocess.pipeline.Convert import Convert
from radprocess import radmc3d
from radprocess import ramses
from radprocess.utils.config import ConfigParams


def _resolve(value, fallback):
    """Return `value` unless it is None, in which case fall back to the config default.

    This is the single convention for method-argument / config interplay:
    every config-backed argument defaults to None in the signature and is
    resolved here, so an explicit call value always wins and an omitted one
    inherits configparams.
    """
    return fallback if value is None else value


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
                         sink_indices=None, gridstyle="octtree", coordsystem="cartesian",
                         min_cells=0):
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
        min_cells : int
            Minimum number of cells in the subbox. Sinks with fewer
            cells are skipped (default 0 = no minimum).

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
            min_cells=min_cells,
        )

    def run_polaris_opacity(
        self,
        dust_mixtures = None,
        mean_molecular_weight = None,
        mass_fraction = None,
        nr_threads = None,
        grid_path = None,
        n_dust_override = None,
        polaris_binary = None,
        cleanup = True,
    ):
        """
        Run POLARIS with 1 photon package to generate dust opacity tables.

        This is Step 4 of the pipeline. It requires that a POLARIS grid
        file already exists (from convert_to_polaris, Step 2).

        Every argument below defaults to None and, when left as None, is
        resolved from the configuration object: the dust setup from
        configparams.dust and the execution settings from configparams.polaris.
        Pass an argument explicitly only to override the config for this call.

        Parameters
        ----------
        dust_mixtures : dict or None
            Dust material definitions. If None, uses configparams.dust.mixtures.
            See DustConfig for the expected structure.
        mean_molecular_weight : float or None
            Gas mean molecular weight. If None, uses configparams.dust.mean_molecular_weight.
        mass_fraction : float or None
            Total dust-to-gas mass fraction. If None, uses configparams.dust.mass_fraction.
        nr_threads : int or None
            Number of OpenMP threads. If None, uses configparams.polaris.nr_threads.
        grid_path : str or Path or None
            Path to the POLARIS grid file. If None, auto-detected.
        n_dust_override : int or None
            Override the number of dust species.
        polaris_binary : str or None
            Name or path of the POLARIS executable. If None, uses
            configparams.polaris.polaris_binary.
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
        dust = self.configparams.dust
        pc = self.configparams.polaris

        # Resolve arguments against the config (explicit call value wins).
        dust_mixtures = _resolve(dust_mixtures, dust.mixtures)
        mean_molecular_weight = _resolve(mean_molecular_weight, dust.mean_molecular_weight)
        mass_fraction = _resolve(mass_fraction, dust.mass_fraction)
        nr_threads = _resolve(nr_threads, pc.nr_threads)
        polaris_binary = _resolve(polaris_binary, pc.polaris_binary)

        if not dust_mixtures:
            raise ValueError(
                "No dust mixtures available: set configparams.dust.mixtures "
                "or pass dust_mixtures=... (see DustConfig for the structure)."
            )

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
            ramses_dir = ramses_dir,
            polaris_dir = polaris_dir,
            grid_path = grid_path,
            dust_mixtures = dust_mixtures,
            mean_molecular_weight = mean_molecular_weight,
            mass_fraction = mass_fraction,
            nr_threads = nr_threads,
            f_acc = f_acc,
            n_dust_override = n_dust_override,
            polaris_binary = polaris_binary,
            cleanup = cleanup,
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
        nphot_scat=None,
        setthreads=None,
        scattering_mode=None,
        scattering_mode_max=None,
        modified_random_walk=None,
        rto_style=None,
        rto_single=None,
        subbox=False,
    ):
        """
        Convert POLARIS opacities to RADMC-3D format and write all
        remaining RADMC-3D input files (Step 5).

        Parameters default to values in configparams.radmc3d when not
        explicitly provided.

        Parameters
        ----------
        polaris_data_dir : str or Path or None
            Path to the POLARIS data/ directory. If None, auto-detected.
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
        subbox : bool
            If True, after writing the shared files to the main radmc3d/
            directory, distribute them (dustkappa_*.inp, dustopac.inp,
            wavelength_micron.inp, radmc3d.inp) to all subbox folders
            in radmc3d/subboxes/. stars.inp is NOT overwritten since
            each subbox has its own star (written by convert_subboxes).

        Returns
        -------
        radmc_dir : Path
            The RADMC-3D directory containing all input files.
        """
        from radprocess.radmc3d.prepare import (
            prepare_radmc3d_inputs, write_radmc3d_control,
        )
        import shutil

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
        if nphot_scat is None:
            nphot_scat = rc.nphot_scat
        if scattering_mode_max is None:
            scattering_mode_max = rc.scattering_mode_max
        if modified_random_walk is None:
            modified_random_walk = rc.modified_random_walk
        if rto_style is None:
            rto_style = rc.rto_style
        if rto_single is None:
            rto_single = rc.rto_single

        if polaris_data_dir is None:
            polaris_data_dir = self.polaris_outputs_dir / "data"
            if not polaris_data_dir.exists():
                raise FileNotFoundError(
                    f"POLARIS data directory not found: {polaris_data_dir}. "
                    "Run run_polaris_opacity() first (Step 4)."
                )

        result = prepare_radmc3d_inputs(
            ramses_dir=ramses_dir,
            radmc_dir=radmc_dir,
            polaris_data_dir=polaris_data_dir,
            n_dust=n_dust,
            f_acc=f_acc,
            wave_min=wave_min,
            wave_max=wave_max,
            n_wavelengths=n_wavelengths,
            nphot=nphot,
            nphot_scat=nphot_scat,
            setthreads=setthreads,
            scattering_mode=scattering_mode,
            scattering_mode_max=scattering_mode_max,
            modified_random_walk=modified_random_walk,
            rto_style=rto_style,
            rto_single=rto_single,
        )

        # Distribute shared files to subbox folders
        if subbox:
            subbox_base = radmc_dir / "subboxes"
            if not subbox_base.exists():
                raise FileNotFoundError(
                    f"Subbox directory not found: {subbox_base}. "
                    "Run convert_subboxes(which_rad='radmc') first."
                )

            # Collect shared files (everything except stars.inp,
            # amr_grid.inp, dust_density.inp which are per-sink)
            shared_files = (
                list(radmc_dir.glob("dustkappa_*.inp"))
                + [radmc_dir / "dustopac.inp",
                   radmc_dir / "wavelength_micron.inp",
                   radmc_dir / "radmc3d.inp"]
            )

            # Read the wavelength grid (same one just written)
            wavelengths = np.logspace(np.log10(wave_min), np.log10(wave_max), n_wavelengths)

            sink_dirs = sorted([
                d for d in subbox_base.iterdir()
                if d.is_dir() and d.name.startswith("sink_")
            ])

            # Read sink catalog for stellar properties
            import csv as csv_mod
            catalog_path = subbox_base / "sink_catalog.csv"
            sink_catalog = {}
            if catalog_path.exists():
                with open(catalog_path, "r") as cf:
                    reader = csv_mod.DictReader(cf)
                    for row in reader:
                        sink_catalog[row["folder"]] = row

            # Read full sink data for stellar property computation
            from radprocess.constants.constants import (
                Ggram, M_sun, L_sun, R_sun, sigma,
            )
            sinks = ramses.read.sink_info(ramses_dir)
            cols = sinks.columns
            m_col = cols.index("M[Msol]")
            accrate_col = cols.index("acc_rate[Msol/y]")
            acclum_col = cols.index("acc_lum[Lsol]")
            intlum_col = cols.index("int_lum[Lsol]")
            teff_col = cols.index("Teff[K]")
            sec_yr = 365.25 * 24 * 3600

            print(f"\nDistributing shared files to {len(sink_dirs)} subbox folders...")
            for sink_dir in sink_dirs:
                # Copy shared files (except radmc3d.inp which is per-sink)
                for src in shared_files:
                    if src.exists():
                        if src.name == "radmc3d.inp":
                            continue  # written per-sink below
                        dst = sink_dir / src.name
                        shutil.copy2(src, dst)

                # Write per-sink stars.inp and radmc3d.inp
                folder_name = sink_dir.name
                if folder_name in sink_catalog:
                    row_idx = int(sink_catalog[folder_name]["row_idx"])

                    sink_mass_g = sinks.data[row_idx, m_col] * M_sun
                    lacc = sinks.data[row_idx, acclum_col] * L_sun
                    lint = sinks.data[row_idx, intlum_col] * L_sun
                    ltot = lint + lacc
                    teff_K = sinks.data[row_idx, teff_col]

                    if lacc > 0 and teff_K > 0:
                        sink_radius_cm = np.sqrt(
                            ltot / (4 * np.pi * sigma * teff_K**4)
                        )
                    else:
                        sink_radius_cm = 1.0 * R_sun

                    if teff_K <= 0:
                        teff_K = 5e3

                    # Read sink offset from box center (in cm for RADMC-3D)
                    offset_file = sink_dir / "sink_offset.txt"
                    if offset_file.exists():
                        offset = np.loadtxt(offset_file)
                        star_x, star_y, star_z = offset[0], offset[1], offset[2]
                    else:
                        star_x, star_y, star_z = 0.0, 0.0, 0.0

                    filepath = sink_dir / "stars.inp"
                    with open(filepath, "w") as sf:
                        sf.write("2\n")
                        sf.write(f"1 {len(wavelengths)}\n")
                        sf.write(
                            f"{sink_radius_cm:.6e} {sink_mass_g:.6e} "
                            f"{star_x:.6e} {star_y:.6e} {star_z:.6e}\n"
                        )
                        for w in wavelengths:
                            sf.write(f"{w:e}\n")
                        sf.write(f"-{teff_K:e}\n")

                    # ---- Per-sink radmc3d.inp with subbox_* parameters ----
                    # Regenerate from config (single source of truth) and add
                    # the subbox_regrid cube centered on the sink at its FOV.
                    # The grid is padded (bigger than the requested FOV) so the
                    # octree is artifact-free. subbox_* are the "long tail"
                    # radmc3d.inp keys, passed via extra_options rather than the
                    # writer's core signature.
                    subbox_opts = {}
                    fov_file = sink_dir / "requested_hw_au.txt"
                    if fov_file.exists():
                        from radprocess.constants.constants import au2cm
                        req_hw_cm = float(np.loadtxt(fov_file)) * au2cm
                        subbox_npix = 128
                        subbox_opts = {
                            "subbox_nx": subbox_npix,
                            "subbox_ny": subbox_npix,
                            "subbox_nz": subbox_npix,
                            "subbox_x0": f"{star_x - req_hw_cm:.6e}",
                            "subbox_x1": f"{star_x + req_hw_cm:.6e}",
                            "subbox_y0": f"{star_y - req_hw_cm:.6e}",
                            "subbox_y1": f"{star_y + req_hw_cm:.6e}",
                            "subbox_z0": f"{star_z - req_hw_cm:.6e}",
                            "subbox_z1": f"{star_z + req_hw_cm:.6e}",
                        }

                    write_radmc3d_control(
                        sink_dir,
                        nphot=nphot,
                        nphot_scat=nphot_scat,
                        setthreads=setthreads,
                        scattering_mode=scattering_mode,
                        scattering_mode_max=scattering_mode_max,
                        modified_random_walk=modified_random_walk,
                        rto_style=rto_style,
                        rto_single=rto_single,
                        extra_options=subbox_opts,
                    )

                print(f"    {sink_dir.name}: OK")

            print(f"Shared files + per-sink stars.inp + radmc3d.inp distributed.\n")

        return result

    def run_radmc3d_mctherm(self, radmc3d_binary="radmc3d", subbox=None):
        """
        Execute ``radmc3d mctherm`` to compute the dust temperature (Step 6).

        All input files must already exist in the target directory
        (from Steps 3 and 5, and distribute via prepare_radmc3d_inputs
        with subbox=True).

        Parameters
        ----------
        radmc3d_binary : str
            Name or path of the RADMC-3D executable.
        subbox : None, str, list of str, or "all"
            Which directory to run mctherm in:
            - None: run on the main radmc3d/ directory (full cloud).
            - "sink_0042": run on a single subbox.
            - ["sink_0003", "sink_0007"]: run on a list of subboxes.
            - "all": run on all subbox folders found in radmc3d/subboxes/.

        Returns
        -------
        result : Path or dict
            If subbox is None or a single string: Path to dust_temperature file.
            If subbox is a list or "all": dict {sink_name: Path} for each.
        """
        from radprocess.radmc3d.run import mctherm

        radmc_dir = self.radmc_outputs_dir

        # No subbox: run on the main directory (original behavior)
        if subbox is None:
            root = self.get_amr_root()
            output_num = int(root.attrs.get("ramses_output_num", 0))
            log_path = radmc_dir / f"radmc3d_mctherm_{output_num:05d}.log"

            return mctherm(
                radmc_dir=radmc_dir,
                log_path=log_path,
                radmc3d_binary=radmc3d_binary,
            )

        # Resolve subbox list
        subbox_base = radmc_dir / "subboxes"
        if not subbox_base.exists():
            raise FileNotFoundError(
                f"Subbox directory not found: {subbox_base}. "
                "Run convert_subboxes(which_rad='radmc') first."
            )

        if subbox == "all":
            sink_names = sorted([
                d.name for d in subbox_base.iterdir()
                if d.is_dir() and d.name.startswith("sink_")
            ])
        elif isinstance(subbox, str):
            sink_names = [subbox]
        elif isinstance(subbox, list):
            sink_names = subbox
        else:
            raise ValueError(
                f"subbox must be None, 'all', a string, or a list of strings, "
                f"got {type(subbox)}"
            )

        if not sink_names:
            raise RuntimeError("No subbox folders found.")

        print(f"\nRunning mctherm on {len(sink_names)} subbox(es)...\n")

        results = {}
        for i, name in enumerate(sink_names):
            sink_dir = subbox_base / name
            if not sink_dir.exists():
                print(f"WARNING: {sink_dir} not found, skipping.")
                continue

            print(f"{'='*60}")
            print(f"  mctherm: {name}  ({i+1}/{len(sink_names)})")
            print(f"{'='*60}")

            log_path = sink_dir / "radmc3d_mctherm.log"

            temp_file = mctherm(
                radmc_dir=sink_dir,
                log_path=log_path,
                radmc3d_binary=radmc3d_binary,
            )
            results[name] = temp_file

        print(f"\nAll mctherm runs completed: {len(results)}/{len(sink_names)} succeeded.\n")

        # If single subbox was passed as string, return just the path
        if isinstance(subbox, str) and subbox != "all":
            return results.get(subbox)

        return results

    def merge_temperature(self, n_dust=None, subbox=None):
        """
        Merge RADMC-3D dust temperatures into the POLARIS grid (Step 7/8).

        For the full-cloud case (subbox=None):
            Reads the POLARIS grid_temp.dat (from Step 4) and the RADMC-3D
            dust_temperature.bdat (from Step 6), replaces the POLARIS dust
            temperatures with the RADMC-3D values, and writes the result as
            grid_temp.radmc3d.dat in the POLARIS output directory.

        For subboxes:
            For each subbox, reads the POLARIS subbox grid
            (ramses_grid_sink_XXXX.dat) and the RADMC-3D temperature
            (dust_temperature.bdat), injects the dust temperatures as
            new parameters (ID=2), and writes grid_temp.radmc3d.dat
            in the POLARIS subbox directory.

        Parameters
        ----------
        n_dust : int or None
            Number of dust species. If None, auto-detected.
        subbox : None, str, list of str, or "all"
            - None: merge on the main directories (full cloud).
            - "sink_0042": merge a single subbox.
            - ["sink_0003", "sink_0007"]: merge a list of subboxes.
            - "all": merge all subbox folders.

        Returns
        -------
        result : Path or dict
            If subbox is None or a single string: Path to merged grid.
            If subbox is a list or "all": dict {sink_name: Path}.
        """
        from radprocess.polaris.merge import merge_radmc3d_temperature, merge_temperature_into_grid

        ramses_dir = self.configparams.dir.ramses_output
        polaris_dir = self.polaris_outputs_dir
        radmc_dir = self.radmc_outputs_dir

        # Full-cloud case
        if subbox is None:
            return merge_radmc3d_temperature(
                polaris_dir=polaris_dir,
                radmc_dir=radmc_dir,
                n_dust=n_dust,
                ramses_dir=ramses_dir,
            )

        # Resolve subbox list
        polaris_subbox_base = polaris_dir / "subboxes"
        radmc_subbox_base = radmc_dir / "subboxes"

        if not polaris_subbox_base.exists():
            raise FileNotFoundError(
                f"POLARIS subbox directory not found: {polaris_subbox_base}. "
                "Run convert_subboxes(which_rad='polaris') first."
            )
        if not radmc_subbox_base.exists():
            raise FileNotFoundError(
                f"RADMC-3D subbox directory not found: {radmc_subbox_base}. "
                "Run convert_subboxes(which_rad='radmc') first."
            )

        if subbox == "all":
            sink_names = sorted([
                d.name for d in polaris_subbox_base.iterdir()
                if d.is_dir() and d.name.startswith("sink_")
            ])
        elif isinstance(subbox, str):
            sink_names = [subbox]
        elif isinstance(subbox, list):
            sink_names = subbox
        else:
            raise ValueError(
                f"subbox must be None, 'all', a string, or a list of strings, "
                f"got {type(subbox)}"
            )

        # Auto-detect n_dust
        if n_dust is None:
            root = self.get_amr_root()
            if "dust_massdensities" in root:
                n_dust = np.array(root["dust_massdensities"]).shape[1]
            elif "dust_massdensity" in root:
                arr = np.array(root["dust_massdensity"])
                n_dust = arr.shape[1] if arr.ndim > 1 else 1
            else:
                n_dust = 1
            print(f"Auto-detected {n_dust} dust species.\n")

        print(f"\nMerging temperatures for {len(sink_names)} subbox(es)...\n")

        results = {}
        for i, name in enumerate(sink_names):
            polaris_sink_dir = polaris_subbox_base / name
            radmc_sink_dir = radmc_subbox_base / name

            if not polaris_sink_dir.exists():
                print(f"WARNING: POLARIS subbox {polaris_sink_dir} not found, skipping.")
                continue
            if not radmc_sink_dir.exists():
                print(f"WARNING: RADMC-3D subbox {radmc_sink_dir} not found, skipping.")
                continue

            # Find POLARIS grid
            polaris_grids = list(polaris_sink_dir.glob("ramses_grid_sink_*.dat"))
            if not polaris_grids:
                print(f"WARNING: No POLARIS grid found in {polaris_sink_dir}, skipping.")
                continue
            polaris_grid_file = polaris_grids[0]

            # Find RADMC-3D temperature
            radmc_temp_file = radmc_sink_dir / "dust_temperature.bdat"
            if not radmc_temp_file.exists():
                radmc_temp_file = radmc_sink_dir / "dust_temperature.dat"
            if not radmc_temp_file.exists():
                print(f"WARNING: No temperature file in {radmc_sink_dir}, skipping.")
                continue

            output_file = polaris_sink_dir / "grid_temp.radmc3d.dat"

            print(f"{'='*60}")
            print(f"  Merging: {name}  ({i+1}/{len(sink_names)})")
            print(f"{'='*60}")

            merged, n_cells = merge_temperature_into_grid(
                polaris_grid_file=polaris_grid_file,
                radmc_temp_file=radmc_temp_file,
                output_file=output_file,
                n_dust=n_dust,
            )
            results[name] = merged

        print(f"\nAll merges completed: {len(results)}/{len(sink_names)} succeeded.\n")

        if isinstance(subbox, str) and subbox != "all":
            return results.get(subbox)

        return results

    def render_images(
        self,
        dust_mixtures = None,
        npix = None,
        distance_pc = None,
        wavelengths_mm = None,
        views = None,
        fov_m = None,
        fov_au = None,
        label = "whole",
        grid_path = None,
        n_dust = None,
        nr_threads = None,
        mean_molecular_weight = None,
        mass_fraction = None,
        polaris_binary = None,
        cleanup_views = True,
        polaris_cmd = None,
        alignment = None,
        peel_off = None,
        acceptance_angle = None,
        nr_photons_scat = None,
        source_star_scat = None,
        subbox = None,
    ):
        """
        Run POLARIS dust emission/scattering imaging.

        The observation setup defaults to configparams.imaging, the dust setup
        to configparams.dust, and the execution settings to configparams.polaris.
        Every config-backed argument defaults to None and is resolved from the
        config; pass one explicitly only to override it for this call.

        Parameters
        ----------
        dust_mixtures : dict or None
            Dust material definitions. If None, uses configparams.dust.mixtures.
        npix : int or None
            Image resolution (npix x npix pixels). If None, configparams.imaging.npix.
        distance_pc : float or None
            Source distance in parsecs. If None, configparams.imaging.distance_pc (required).
        wavelengths_mm : list of float or None
            Wavelengths to image in millimetres. If None,
            configparams.imaging.wavelengths_mm (required, non-empty).
        views : list of str or None
            Viewing angles (subset of ["xy", "xz", "yz"]). If None,
            configparams.imaging.views.
        fov_m : float or None
            Field of view in metres (full width). Overrides fov_au if set.
        fov_au : float or None
            Field of view in AU (full width). If None, configparams.imaging.fov_au.
            For subboxes, use the same value as the RADMC-3D regrid box_au.
        label : str
            Output subdirectory label (per-call, not config-backed).
        grid_path : str or Path or None
            Path to the merged grid. If None, auto-detected.
        n_dust : int or None
            Number of dust species. If None, auto-detected.
        nr_threads : int or None
            OpenMP threads. If None, configparams.polaris.nr_threads.
        mean_molecular_weight : float or None
            Gas mu. If None, configparams.dust.mean_molecular_weight.
        mass_fraction : float or None
            Total dust-to-gas mass fraction. If None, configparams.dust.mass_fraction.
        polaris_binary : str or None
            POLARIS executable. If None, configparams.polaris.polaris_binary.
        cleanup_views : bool
            Remove previous image outputs (per-call).
        polaris_cmd : str or None
            POLARIS command. If None, configparams.imaging.polaris_cmd.
        alignment : str or None
            Grain alignment mechanism. If None, configparams.imaging.alignment.
        peel_off : bool or None
            Use peel-off technique. If None, configparams.imaging.peel_off.
        acceptance_angle : float or None
            Acceptance angle for scattered light. If None, configparams.imaging.acceptance_angle.
        nr_photons_scat : int or None
            Photon packages for scattering MC. If None, configparams.imaging.nr_photons_scat.
        source_star_scat : list of dict or None
            Stellar sources for scattering (per-call).
        subbox : None, str, list of str, or "all"
            - None: render full cloud.
            - "sink_0042": single subbox.
            - ["sink_0003", "sink_0007"]: list of subboxes.
            - "all": all subboxes with merged grids.

        Returns
        -------
        result : dict
            If subbox is None: {view_name: Path}.
            If subbox is set: {sink_name: {view_name: Path}}.
        """
        from radprocess.polaris.imaging import render_images

        polaris_dir = self.polaris_outputs_dir
        img = self.configparams.imaging
        dust = self.configparams.dust
        pc = self.configparams.polaris

        # Resolve arguments against the config (explicit call value wins).
        dust_mixtures = _resolve(dust_mixtures, dust.mixtures)
        npix = _resolve(npix, img.npix)
        distance_pc = _resolve(distance_pc, img.distance_pc)
        wavelengths_mm = _resolve(wavelengths_mm, img.wavelengths_mm)
        views = _resolve(views, img.views)
        fov_au = _resolve(fov_au, img.fov_au)
        polaris_cmd = _resolve(polaris_cmd, img.polaris_cmd)
        alignment = _resolve(alignment, img.alignment)
        peel_off = _resolve(peel_off, img.peel_off)
        acceptance_angle = _resolve(acceptance_angle, img.acceptance_angle)
        nr_photons_scat = _resolve(nr_photons_scat, img.nr_photons_scat)
        mean_molecular_weight = _resolve(mean_molecular_weight, dust.mean_molecular_weight)
        mass_fraction = _resolve(mass_fraction, dust.mass_fraction)
        nr_threads = _resolve(nr_threads, pc.nr_threads)
        polaris_binary = _resolve(polaris_binary, pc.polaris_binary)

        # Required fields (no safe default): fail early with a clear message.
        if not dust_mixtures:
            raise ValueError(
                "No dust mixtures available: set configparams.dust.mixtures "
                "or pass dust_mixtures=... (see DustConfig)."
            )
        if distance_pc is None:
            raise ValueError(
                "distance_pc is required: set configparams.imaging.distance_pc "
                "or pass distance_pc=..."
            )
        if not wavelengths_mm:
            raise ValueError(
                "wavelengths_mm is required: set configparams.imaging.wavelengths_mm "
                "or pass wavelengths_mm=..."
            )

        # Honour a config/argument fov_au in the full-cloud path too (the
        # subbox path handles fov_au/fov_m itself further below).
        if fov_m is None and fov_au is not None:
            from radprocess.constants.constants import au2m as _au2m
            fov_m = fov_au * _au2m

        if n_dust is None:
            root = self.get_amr_root()
            if "dust_massdensities" in root:
                n_dust = np.array(root["dust_massdensities"]).shape[1]
            elif "dust_massdensity" in root:
                arr = np.array(root["dust_massdensity"])
                n_dust = arr.shape[1] if arr.ndim > 1 else 1
            else:
                n_dust = 1

        root = self.get_amr_root()
        output_num = int(root.attrs.get("ramses_output_num", 0))

        # ---- Full-cloud case ----
        if subbox is None:
            image_output_dir = polaris_dir / "images"

            if grid_path is None:
                grid_path = polaris_dir / "grid_temp.radmc3d.dat"
                if not grid_path.exists():
                    raise FileNotFoundError(
                        f"Merged grid not found: {grid_path}. "
                        "Run merge_temperature() first."
                    )

            return render_images(
                polaris_dir = polaris_dir,
                image_output_dir = image_output_dir,
                grid_path = grid_path,
                dust_mixtures = dust_mixtures,
                n_dust = n_dust,
                mean_molecular_weight = mean_molecular_weight,
                mass_fraction = mass_fraction,
                npix = npix,
                distance_pc = distance_pc,
                wavelengths_mm = wavelengths_mm,
                views = views,
                nr_threads = nr_threads,
                fov_m = fov_m,
                output_num = output_num,
                polaris_binary = polaris_binary,
                label = label,
                cleanup_views = cleanup_views,
                polaris_cmd = polaris_cmd,
                alignment = alignment,
                peel_off = peel_off,
                acceptance_angle = acceptance_angle,
                nr_photons_scat = nr_photons_scat,
                source_star_scat = source_star_scat,
            )

        # ---- Subbox case ----
        polaris_subbox_base = polaris_dir / "subboxes"
        if not polaris_subbox_base.exists():
            raise FileNotFoundError(
                f"POLARIS subbox directory not found: {polaris_subbox_base}. "
                "Run convert_subboxes(which_rad='polaris') first."
            )

        if subbox == "all":
            sink_names = sorted([
                d.name for d in polaris_subbox_base.iterdir()
                if d.is_dir() and d.name.startswith("sink_")
                and (d / "grid_temp.radmc3d.dat").exists()
            ])
        elif isinstance(subbox, str):
            sink_names = [subbox]
        elif isinstance(subbox, list):
            sink_names = subbox
        else:
            raise ValueError(
                f"subbox must be None, 'all', a string, or a list, "
                f"got {type(subbox)}"
            )

        if not sink_names:
            raise RuntimeError("No subbox folders with merged grids found.")

        print(f"\nRendering images for {len(sink_names)} subbox(es)...\n")

        all_results = {}
        for i, name in enumerate(sink_names):
            polaris_sink_dir = polaris_subbox_base / name

            if not polaris_sink_dir.exists():
                print(f"WARNING: {polaris_sink_dir} not found, skipping.")
                continue

            sink_grid = polaris_sink_dir / "grid_temp.radmc3d.dat"
            if not sink_grid.exists():
                print(f"WARNING: No merged grid in {polaris_sink_dir}, skipping.")
                continue

            print(f"\n{'='*60}")
            print(f"  Rendering: {name}  ({i+1}/{len(sink_names)})")
            print(f"{'='*60}")

            image_output_dir = polaris_dir / "images" / "subboxes" / name

            # Read sink offset (in metres) for detector centering
            sink_offset_m = None
            offset_file = polaris_sink_dir / "sink_offset.txt"
            if offset_file.exists():
                sink_offset_m = np.loadtxt(offset_file).tolist()

            # Determine the FOV for this subbox
            from radprocess.constants.constants import au2m as _au2m
            if fov_au is not None:
                subbox_fov_m = fov_au * _au2m
            elif fov_m is not None:
                subbox_fov_m = fov_m
            else:
                # Fall back to the requested half-width from extraction
                fov_file = polaris_sink_dir / "requested_hw_au.txt"
                if fov_file.exists():
                    req_hw_au = float(np.loadtxt(fov_file))
                    subbox_fov_m = req_hw_au * _au2m * 2.0
                else:
                    subbox_fov_m = None

            if subbox_fov_m is not None:
                print(f"    FOV: {subbox_fov_m / _au2m:.0f} AU, centered on sink")

            image_dirs = render_images(
                polaris_dir = polaris_sink_dir,
                image_output_dir = image_output_dir,
                grid_path = sink_grid,
                dust_mixtures = dust_mixtures,
                n_dust = n_dust,
                mean_molecular_weight = mean_molecular_weight,
                mass_fraction = mass_fraction,
                npix = npix,
                distance_pc = distance_pc,
                wavelengths_mm = wavelengths_mm,
                views = views,
                nr_threads = nr_threads,
                fov_m = subbox_fov_m,
                output_num = output_num,
                polaris_binary = polaris_binary,
                label = label,
                cleanup_views = cleanup_views,
                polaris_cmd = polaris_cmd,
                alignment = alignment,
                peel_off = peel_off,
                acceptance_angle = acceptance_angle,
                nr_photons_scat = nr_photons_scat,
                source_star_scat = source_star_scat,
                detector_shift_m = sink_offset_m,
            )

            all_results[name] = image_dirs

        print(f"\nAll rendering completed: {len(all_results)}/{len(sink_names)} succeeded.\n")

        if isinstance(subbox, str) and subbox != "all":
            return all_results.get(subbox)

        return all_results