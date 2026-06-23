from dataclasses import dataclass, fields, is_dataclass, field
#from pydantic import BaseModel, Field
import html


def validate_type(obj, field_name, value):
    """Validate the type of a field dynamically based on the dataclass annotation."""
    flds = fields(obj)
    f = next(ff for ff in flds if ff.name == field_name)
    expected_type = f.type

    if not isinstance(value, expected_type):
        raise TypeError(
            f"Invalid type for '{field_name}': "
            f"expected {expected_type.__name__}, got {type(value).__name__}"
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
        typ = f.type.__name__ if hasattr(f.type, "__name__") else str(f.type)

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

def _repr_html_(self):
    body = (
        _html_table(self.inout)
        + "<br/>"
        + _html_table(self.amrsource)
        + "<br/>"
        + _html_table(self.sim)
    )

    return f"""
    <details style="margin:8px 0; padding:6px; border:1px solid #ccc; border-radius:6px;">
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


def __repr__(self):
    return (
        "PARAMS\n-------\n"
        + fancy_repr(self.inout)
        + "\n\n"
        + fancy_repr(self.sim)
    )

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
    max_type = max(len(f.type.__name__) if hasattr(f.type, "__name__") else len(str(f.type))
                   for f in flds)
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
        typ = (
            f.type.__name__.ljust(max_type)
            if hasattr(f.type, "__name__")
            else str(f.type).ljust(max_type)
        )
        desc = f.metadata.get("desc", "").ljust(max_desc)
        value = fmt_value(getattr(obj, f.name)).ljust(VALUE_COL_WIDTH)
        lines.append(f"{name}  {typ}  {value}  {desc}")

    return "\n".join(lines)


class StrictDataclass:
    def __setattr__(self, name, value):
        if name in {f.name for f in fields(self)}:
            expected_type = next(f for f in fields(self) if f.name == name).type
            if not isinstance(value, expected_type):
                raise TypeError(
                    f"❌ Error: Field '{name}' expects type {expected_type.__name__}, "
                    f"got {type(value).__name__}"
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
    dust_size_min: float = field(default=5e-9,
        metadata={'desc': r'[m] Minimum grain radius.'})
    dust_size_max: float = field(default=2.5e-7,
        metadata={'desc': r'[m] Maximum grain radius.'})
    dust_size_powerlaw: float = field(default=-3.5,
        metadata={'desc': r'Power-law exponent for the grain size distribution (e.g. -3.5 for MRN).'})
    mean_molecular_weight: float = field(default=2.37,
        metadata={'desc': r'Mean molecular weight (mu) of the gas.'})
    mass_fraction: float = field(default = 1,
        metadata={'desc': r'Dust-to-gas mass fraction.'})
    nr_threads: int = field(default=8,
        metadata={'desc': r'Number of OpenMP threads for POLARIS.'})
    polaris_binary: str = field(default='polaris',
        metadata={'desc': r'Name or path of the POLARIS executable.'})

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
    radmc3d: Radmc3dConfig = field(default_factory=Radmc3dConfig)
    polaris: PolarisConfig = field(default_factory=PolarisConfig)
    nb_dust: int = 0

    def __repr__(self):
        return (
            "PARAMS\n-------\n"
            + fancy_repr(self.dir)
            + "\n\n"
            + fancy_repr(self.amrsource)
            + "\n\n"
            + fancy_repr(self.sim)
            + "\n\n"
            + fancy_repr(self.radmc3d)
            + "\n\n"
            + fancy_repr(self.polaris)
        )

    def _repr_html_(self):
        body = (
            _html_table(self.dir)
            + "<br/>"
            + _html_table(self.amrsource)
            + "<br/>"
            + _html_table(self.sim)
            + "<br/>"
            + _html_table(self.radmc3d)
            + "<br/>"
            + _html_table(self.polaris)
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
