from dataclasses import dataclass, fields, is_dataclass, field
from typing import Optional, Union, get_origin, get_args
#from pydantic import BaseModel, Field
import html


# =======================================================
# Type helpers (support plain types, Optional/Union, and
# bare containers such as list / dict)
# =======================================================
def _type_name(expected):
    """Human-readable name for a type annotation (for repr / error messages)."""
    origin = get_origin(expected)
    if origin is Union:
        args = get_args(expected)
        non_none = [a for a in args if a is not type(None)]
        names = " | ".join(_type_name(a) for a in non_none)
        if type(None) in args:
            return f"Optional[{names}]" if len(non_none) == 1 else f"{names} | None"
        return names
    if hasattr(expected, "__name__"):
        return expected.__name__
    return str(expected).replace("typing.", "")


def _check_type(value, expected):
    """Return True if `value` matches the annotation `expected`.

    Handles:
      * plain types / dataclasses  -> isinstance check
      * Optional[X] / Union[...]   -> matches any member (incl. NoneType)
      * parameterized generics     -> checks the container origin only
        (e.g. list[float] -> isinstance(value, list))
    """
    origin = get_origin(expected)
    if origin is None:
        return isinstance(value, expected)
    if origin is Union:
        return any(_check_type(value, arg) for arg in get_args(expected))
    # list[...], dict[...], tuple[...], etc.: validate the container type only
    return isinstance(value, origin)


def validate_type(obj, field_name, value):
    """Validate the type of a field dynamically based on the dataclass annotation."""
    flds = fields(obj)
    f = next(ff for ff in flds if ff.name == field_name)
    expected_type = f.type

    if not _check_type(value, expected_type):
        raise TypeError(
            f"Invalid type for '{field_name}': "
            f"expected {_type_name(expected_type)}, got {type(value).__name__}"
        )


# =======================================================
# Utility: scientific number formatting
# =======================================================
def fmt_value(v, sci_threshold=1e3):
    # Booleans first
    if isinstance(v, bool):
        return "True" if v else "False"

    # Numbers
    if isinstance(v, (float, int)):
        if v != 0 and (abs(v) >= sci_threshold or abs(v) < 1e10):
            return f"{v:.3e}"

    return repr(v)

def html_value(v):
    """Return HTML representation with boolean badges."""
    if isinstance(v, bool):
        if v:
            return (
                '<span style="background:#c6f6c6; color:#005000; '
                'padding:2px 6px; border-radius:6px; font-weight:bold;'
                'border:1px solid #8ad88a;">True</span>'
            )
        else:
            return (
                '<span style="background:#f8c8c8; color:#7a0000; '
                'padding:2px 6px; border-radius:6px; font-weight:bold;'
                'border:1px solid #e08a8a;">False</span>'
            )

    # Default: code styling for numbers & strings
    return f"<code>{fmt_value(v)}</code>"

# =======================================================
# HTML REPR FOR JUPYTER (collapsible panel)
# =======================================================
def _html_table(obj):
    """Return only the table HTML for a dataclass (no <details>)."""
    cls_name = obj.__class__.__name__
    flds = fields(obj)

    rows = ""
    for f in flds:
        desc_raw = f.metadata.get("desc", "")
        desc = html.escape(desc_raw).replace("&dollar;", "$")

        val = getattr(obj, f.name)
        typ = _type_name(f.type)

        rows += f"""
        <tr>
            <td style="text-align:left;"><b>{f.name}</b></td>
            <td style="text-align:left;">{typ}</td>
            <td style="text-align:left;">{html_value(val)}</td>
            <td style="text-align:left; white-space:normal;">{desc}</td>
        </tr>
        """

    return f"""
    <h4 style="margin-bottom:3px; margin-top:18px; text-align:left;">{cls_name}</h4>

    <table style="border-collapse: collapse; margin-top:5px; text-align:left;">
      <thead>
        <tr style="text-align:left; border-bottom:1px solid #888;">
          <th style="padding-right:20px; text-align:left;">Field</th>
          <th style="padding-right:20px; text-align:left;">Type</th>
          <th style="padding-right:20px; text-align:left;">Value</th>
          <th style="text-align:left;">Description</th>
        </tr>
      </thead>
      <tbody>
        {rows}
      </tbody>
    </table>
    """



def _html_repr(obj):
    """Collapsible wrapper for single dataclasses (collapsed by default)."""
    return f"""
    <details style="margin:8px 0; padding:6px; border:1px solid #ccc; border-radius:6px;" >
      <summary style="font-size:16px; font-weight:bold; cursor:pointer;">
          {obj.__class__.__name__}
      </summary>
      <div style="padding:10px;">
        {_html_table(obj)}
      </div>
    </details>

    <script>
    if (window.MathJax) {{
        MathJax.typesetPromise();
    }}
    </script>
    """


# =======================================================
# TEXT FALLBACK (for terminal printing)
# =======================================================
def fancy_repr(obj) -> str:
    """Non-HTML pretty repr for terminal."""
    if not is_dataclass(obj):
        return repr(obj)

    cls_name = obj.__class__.__name__
    flds = fields(obj)

    max_name = max(len(f.name) for f in flds)
    max_type = max(len(_type_name(f.type)) for f in flds)
    max_desc = max(len(f.metadata.get("desc", "")) for f in flds)

    VALUE_COL_WIDTH = 14

    lines = [f"{cls_name}:"]
    lines.append(
        f"{'Field'.ljust(max_name)}  {'Type'.ljust(max_type)}  "
        f"{'Value'.ljust(VALUE_COL_WIDTH)}  {'Description'.ljust(max_desc)}"
    )
    lines.append(
        f"{'-'*max_name}  {'-'*max_type}  {'-'*VALUE_COL_WIDTH}  {'-'*max_desc}"
    )

    for f in flds:
        name = f.name.ljust(max_name)
        typ = _type_name(f.type).ljust(max_type)
        desc = f.metadata.get("desc", "").ljust(max_desc)
        value = fmt_value(getattr(obj, f.name)).ljust(VALUE_COL_WIDTH)
        lines.append(f"{name}  {typ}  {value}  {desc}")

    return "\n".join(lines)


class StrictDataclass:
    def __setattr__(self, name, value):
        known_fields = {f.name for f in fields(self)}
        if name in known_fields:
            expected_type = next(f for f in fields(self) if f.name == name).type
            if not _check_type(value, expected_type):
                raise TypeError(
                    f"Error: Field '{name}' expects type {_type_name(expected_type)}, "
                    f"got {type(value).__name__}"
                )
        else:
            raise AttributeError(
                f"Error: '{type(self).__name__}' has no parameter '{name}'. "
                f"Valid parameters are: {', '.join(sorted(known_fields))}"
            )
        super().__setattr__(name, value)



# =======================================================
# PARAMETER DATACLASSES
# =======================================================
@dataclass
class Directories(StrictDataclass):
    ramses_output: str = field(default='ramses_outputs/', 
        metadata={'desc': r'The RAMSES simulation output directory path.'})
    pipeline_output: str = field(default='pipeline_outputs/', 
        metadata={'desc': r'The main pipeline output directory path.'})
    
    def __repr__(self):  # terminal
        return fancy_repr(self)

    def _repr_html_(self):  # Jupyter
        return _html_repr(self)


@dataclass
class AmrSource(StrictDataclass):
    rho: bool = field(default=True, 
        metadata={"ramses_name": "density", "type": "scalar", 'desc': r'Include the gas density field (always True)'})
    fluiddens: bool = field(default=False, 
        metadata={"ramses_name": "fluid_density", "type": "scalar", 'desc': r'Include the multiple fluid densities (if they exist)'})
    dustratios: bool = field(default=False, 
        metadata={"ramses_name": "dust_ratio", "type": "scalar", 'desc': r'Include the multiple dust ratios (if they exist)'})
    vel: bool = field(default=False, 
        metadata={"ramses_name": "velocity", "type": "vector", 'desc': r'Include the velocity field'})
    fluidvel: bool = field(default=False, 
        metadata={"ramses_name": "fluid_v", "type": "vector", 'desc': r'Include the fluid velocity field (if they exist)'})
    bl: bool = field(default=False, 
        metadata={"ramses_name": "B_left", "type": "vector", 'desc': r'Include the left magnetic field (B) components'})
    br: bool = field(default=False, 
        metadata={"ramses_name": "B_right", "type": "vector",'desc': r'Include the right magnetic field (B) components'})
    p: bool = field(default=False, 
        metadata={"ramses_name": "thermal_pressure", "type": "scalar", 'desc': r'Include the pressure field'})
    xi: bool = field(default=False, 
        metadata={"ramses_name": "radiative_energy",  "type": "scalar", 'desc': r'Include the ionization fraction field'})
    phi: bool = field(default=False, 
        metadata={"ramses_name": "passive_scalar",  "type": "scalar", 'desc': r'Include the gravitational potential field'})
    g: bool = field(default=False, 
        metadata={"ramses_name": "passive_scalar",  "type": "scalar", 'desc': r'Include the gravitational acceleration field'})
    temp: bool = field(default=False, 
        metadata={"ramses_name": "temperature", "type": "scalar", 'desc': r'Include the gas temperature field'})

    def __repr__(self):  # terminal
        return fancy_repr(self)

    def _repr_html_(self):  # Jupyter
        return _html_repr(self)


@dataclass
class Sim(StrictDataclass):
    size_hole_au: float = field(default=0., 
        metadata={'desc': r'[AU] Size of the central hole to exctrude around the sinks.'})
    dtogas: float = field(default=0.01, 
        metadata={'desc': r'Dust-to-gas mass ratio. Used only if this is not multi-grain mode or if no dust ratios exist in RAMSES output.'})
    facc: float = field(default=0.1, 
        metadata={'desc': r' Accretion fraction converted into radiation. Change at your own risk!'})
    use_ramses_T: bool = field(default=True, 
        metadata={'desc': r' Use RAMSES stellar temperatures (if available in sink info) as stellar inputs for the RT simulations.'})
    use_ramses_acc_rate: bool = field(default=True, 
        metadata={'desc': r' Use RAMSES accretion rates (if available in sink info) to derive stellar radii for the RT simulations.'})
    use_multi_grain: bool = field(default=True, 
        metadata={'desc': r'Do RT with multiple bins if True (if RAMSES output has multiple bins). If False, then dust density is computed using dtogas*rho.'})



    def __repr__(self):  # terminal
        return fancy_repr(self)

    def _repr_html_(self):  # Jupyter
        return _html_repr(self)


# ---------------- Dust properties (shared by opacity + rendering) ----------------
@dataclass
class DustConfig(StrictDataclass):
    """Dust composition and gas medium, shared by run_polaris_opacity and render_images.

    This is the single source of truth for the dust setup, replacing the
    ``dust_mixtures`` dict that used to be passed separately to both the opacity
    and the imaging steps (which risked the two disagreeing).

    Expected ``mixtures`` structure::

        mixtures = {
            "<mixture_id>": {
                "<component_name>": {
                    "path": str,          # optical data file (.nk / .cs)
                    "distribution": str,  # size distribution keyword, e.g. "plaw"
                    "fraction": float,    # species mass fraction WITHIN the mixture
                                          #   (components of a mixture sum to 1)
                    "density": float,     # grain bulk density [kg/m^3]
                    "amin": float,        # [m] minimum grain radius
                    "amax": float,        # [m] maximum grain radius
                    "index": int | list,  # POLARIS size-distribution parameter(s)
                },
            },
        }

    Note the distinction:
      * ``fraction`` (per component)  -> composition, splits dust among species.
      * ``mass_fraction`` (below)     -> total dust-to-gas mass ratio (POLARIS <mass_fraction>).
    """
    mixtures: dict = field(default_factory=dict,
        metadata={'desc': r'Dust mixtures {mixture_id: {component: {path, distribution, fraction, density, amin, amax, index}}}. Required; see class docstring.'})
    mean_molecular_weight: float = field(default=2.37,
        metadata={'desc': r'Mean molecular weight (mu) of the gas.'})
    # NOTE: default was inconsistent in the old code (1.0 in run_opacity vs 0.01 in
    # the wrapper). Consolidated to the ISM-conventional dust-to-gas ratio 0.01.
    # Change here if your setup expects a different value.
    mass_fraction: float = field(default=0.01,
        metadata={'desc': r'Total dust-to-gas mass fraction (written as POLARIS <mass_fraction>).'})

    def validate(self):
        if not self.mixtures:
            raise ValueError(
                "DustConfig.mixtures is empty: define at least one dust mixture "
                "before running the POLARIS opacity or rendering steps."
            )

    def __repr__(self):
        return fancy_repr(self)

    def _repr_html_(self):
        return _html_repr(self)


@dataclass
class Radmc3dConfig(StrictDataclass):
    nphot: int = field(default=1_000_000,
        metadata={'desc': r'Number of photon packages for mctherm.'})
    nphot_scat: int = field(default=1_000_000,
        metadata={'desc': r'Number of photon packages for scattering Monte Carlo.'})
    setthreads: int = field(default=8,
        metadata={'desc': r'Number of OpenMP threads for RADMC-3D.'})
    scattering_mode: int = field(default=1,
        metadata={'desc': r'Scattering mode (1=isotropic, 2=anisotropic with HG, 5=full).'})
    scattering_mode_max: int = field(default=1,
        metadata={'desc': r'Maximum scattering mode.'})
    modified_random_walk: int = field(default=1,
        metadata={'desc': r'Enable modified random walk (1=yes, 0=no). Accelerates optically thick regions.'})
    rto_style: int = field(default=3,
        metadata={'desc': r'Output style for dust_temperature (1=ascii, 3=binary).'})
    rto_single: int = field(default=1,
        metadata={'desc': r'Single precision output (1=yes, 0=no).'})
    wave_min: float = field(default=0.27,
        metadata={'desc': r'[micron] Minimum wavelength for the wavelength grid.'})
    wave_max: float = field(default=3000.0,
        metadata={'desc': r'[micron] Maximum wavelength for the wavelength grid.'})
    n_wavelengths: int = field(default=200,
        metadata={'desc': r'Number of log-spaced wavelength points.'})

    def __repr__(self):
        return fancy_repr(self)

    def _repr_html_(self):
        return _html_repr(self)


@dataclass
class PolarisConfig(StrictDataclass):
    """POLARIS execution settings only.

    Dust composition (mixtures, mu, mass_fraction) now lives in DustConfig, and
    the old dust_size_min/max/powerlaw fields have been removed: grain sizes come
    from each component's amin/amax/distribution in DustConfig.mixtures.
    """
    nr_threads: int = field(default=8,
        metadata={'desc': r'Number of OpenMP threads for POLARIS.'})
    polaris_binary: str = field(default='polaris',
        metadata={'desc': r'Name or path of the POLARIS executable.'})

    def __repr__(self):
        return fancy_repr(self)

    def _repr_html_(self):
        return _html_repr(self)


# ---------------- Imaging / synthetic observation ----------------
@dataclass
class ImagingConfig(StrictDataclass):
    """Definition of the synthetic observation produced by render_images.

    distance_pc and wavelengths_mm have no safe default and are required; they
    are left unset (None / empty) and checked by validate() before a run rather
    than given an arbitrary value that would silently rescale the output.
    """
    npix: int = field(default=256,
        metadata={'desc': r'Image resolution (npix x npix pixels).'})
    distance_pc: Optional[float] = field(default=None,
        metadata={'desc': r'Source distance in parsecs. REQUIRED.'})
    wavelengths_mm: list = field(default_factory=list,
        metadata={'desc': r'Wavelengths to image, in millimetres. REQUIRED (non-empty list).'})
    views: list = field(default_factory=lambda: ["xy", "xz", "yz"],
        metadata={'desc': r'Default viewing angles to render. Built-in: xy, xz, yz; extend via custom_views.'})
    custom_views: dict = field(default_factory=dict,
        metadata={'desc': r'User-defined views merged over the built-ins: {name: {plane_id, axis1, axis2, theta, phi}}.'})
    fov_au: Optional[float] = field(default=None,
        metadata={'desc': r'[AU] Field of view (full width). None = full grid extent.'})
    polaris_cmd: str = field(default="CMD_DUST_EMISSION",
        metadata={'desc': r'POLARIS command: CMD_DUST_EMISSION, CMD_DUST_SCATTERING, or both.'})
    alignment: str = field(default="ALIG_PA",
        metadata={'desc': r'Grain alignment: ALIG_PA, ALIG_IDG, ALIG_RAT, ALIG_INTERNAL.'})
    peel_off: bool = field(default=True,
        metadata={'desc': r'Use the peel-off technique for scattering.'})
    acceptance_angle: Optional[float] = field(default=None,
        metadata={'desc': r'[deg] Acceptance angle for scattered light. None = POLARIS default.'})
    nr_photons_scat: Optional[int] = field(default=None,
        metadata={'desc': r'Number of photon packages for the scattering Monte Carlo. None = no scattering source.'})
    scat_source_radius_rsun: float = field(default=1.0,
        metadata={'desc': r'[Rsun] Radius of the default scattering source star (used when nr_photons_scat is set but no explicit source is given).'})
    scat_source_temp_k: float = field(default=5000.0,
        metadata={'desc': r'[K] Temperature of the default scattering source star.'})

    def validate(self):
        if self.distance_pc is None:
            raise ValueError("ImagingConfig.distance_pc is required (source distance in pc).")
        if not self.wavelengths_mm:
            raise ValueError("ImagingConfig.wavelengths_mm is required (list of wavelengths in mm).")

    def __repr__(self):
        return fancy_repr(self)

    def _repr_html_(self):
        return _html_repr(self)


# ---------------- Structure Params ----------------
@dataclass
class ConfigParams(StrictDataclass):
    dir: Directories = field(default_factory=Directories)
    amrsource: AmrSource = field(default_factory=AmrSource)
    sim: Sim = field(default_factory=Sim)
    dust: DustConfig = field(default_factory=DustConfig)
    radmc3d: Radmc3dConfig = field(default_factory=Radmc3dConfig)
    polaris: PolarisConfig = field(default_factory=PolarisConfig)
    imaging: ImagingConfig = field(default_factory=ImagingConfig)
    nb_dust: int = 0

    def validate(self):
        """Check cross-cutting required fields before a run."""
        self.dust.validate()
        self.imaging.validate()

    def __repr__(self):
        return (
            "PARAMS\n-------\n"
            + fancy_repr(self.dir)
            + "\n\n"
            + fancy_repr(self.amrsource)
            + "\n\n"
            + fancy_repr(self.sim)
            + "\n\n"
            + fancy_repr(self.dust)
            + "\n\n"
            + fancy_repr(self.radmc3d)
            + "\n\n"
            + fancy_repr(self.polaris)
            + "\n\n"
            + fancy_repr(self.imaging)
        )

    def _repr_html_(self):
        body = (
            _html_table(self.dir)
            + "<br/>"
            + _html_table(self.amrsource)
            + "<br/>"
            + _html_table(self.sim)
            + "<br/>"
            + _html_table(self.dust)
            + "<br/>"
            + _html_table(self.radmc3d)
            + "<br/>"
            + _html_table(self.polaris)
            + "<br/>"
            + _html_table(self.imaging)
        )

        return f"""
        <details style="margin:8px 0; padding:6px; border:1px solid #ccc; border-radius:6px;" open>
          <summary style="font-size:16px; font-weight:bold; cursor:pointer;">Configuration Parameters</summary>
          <div style="padding:10px;">
            {body}
          </div>
        </details>

        <script>
        if (window.MathJax) {{
            MathJax.typesetPromise();
        }}
        </script>
        """