from dataclasses import dataclass, fields, is_dataclass, field
#from pydantic import BaseModel, Field
import html
from typing import Literal, Any
import ipywidgets as widgets
from IPython.display import display 


def interactive_config(dataclass_obj):
    widgets_dict = {}

    for f in fields(dataclass_obj):
        val = getattr(dataclass_obj, f.name)
        typ = f.type

        # Boolean -> Checkbox
        if typ is bool:
            w = widgets.Checkbox(value=val, description=f.name)

        # Integer -> IntText
        elif typ is int:
            w = widgets.IntText(value=val, description=f.name)

        # Float -> FloatText
        elif typ is float:
            w = widgets.FloatText(value=val, description=f.name)

        # String -> Text input
        elif typ is str:
            w = widgets.Text(value=val, description=f.name)

        # Fallback
        else:
            w = widgets.Text(value=str(val), description=f"{f.name} (untyped)")

        def make_handler(field_name, widget):
            def handler(change):
                setattr(dataclass_obj, field_name, change["new"])
            return handler

        w.observe(make_handler(f.name, w), "value")
        widgets_dict[f.name] = w

    box = widgets.VBox(list(widgets_dict.values()))
    display(box)


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
        + _html_table(self.pymsesrc)
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
class RamsesOutput(StrictDataclass):
    ramses_output_dir: str = field(default='ramses_outputs/', 
        metadata={'desc': r'The RAMSES output directory path.'})
    
    def __repr__(self):  # terminal
        return fancy_repr(self)

    def _repr_html_(self):  # Jupyter
        return _html_repr(self)


@dataclass
class RadiativeOutput(StrictDataclass):
    polaris_output_dir: str = field(default='polaris_outputs/', 
        metadata={'desc': r'The POLARIS output directory path.'})
    radmc_output_dir: str = field(default='radmc3d_outputs/', 
        metadata={'desc': r'The RADMC3D output directory path.'})
    
    def __repr__(self):  # terminal
        return fancy_repr(self)

    def _repr_html_(self):  # Jupyter
        return _html_repr(self)


@dataclass
class Pymsesrc(StrictDataclass):
    rho: bool = field(default=True, 
        metadata={'desc': r'Include the gas density field (always True)'})
    dustratios: bool = field(default=True, 
        metadata={'desc': r'RAMSES simulation has multiple dust species'})
    vel: bool = field(default=True, 
        metadata={'desc': r'Include the velocity field'})
    bl: bool = field(default=False, 
        metadata={'desc': r'Include the magnetic field (B) components'})
    br: bool = field(default=False, 
        metadata={'desc': r'Include the radial magnetic field (B) components'})
    p: bool = field(default=False, 
        metadata={'desc': r'Include the pressure field'})
    xi: bool = field(default=False, 
        metadata={'desc': r'Include the ionization fraction field'})
    phi: bool = field(default=False, 
        metadata={'desc': r'Include the gravitational potential field'})
    g: bool = field(default=False, 
        metadata={'desc': r'Include the gravitational acceleration field'})

    def __repr__(self):  # terminal
        return fancy_repr(self)

    def _repr_html_(self):  # Jupyter
        return _html_repr(self)


@dataclass
class Sim(StrictDataclass):
    size_hole_au: float = field(default=4., 
        metadata={'desc': r'[AU] Size of the central hole to exclude around the star.'})
    dtogas: float = field(default=0.01, 
        metadata={'desc': r'Dust-to-gas mass ratio. Uused only if multi_grain_mode==False or if no dust ratios exist in RAMSES output.'})
    facc: float = field(default=0.1, 
        metadata={'desc': r' Accretion fraction converted into radiation. Change at your own risk!'})
    use_ramses_T: bool = field(default=True, 
        metadata={'desc': r' Use RAMSES stellar temperature field(s) directly (if available) to create the star(s) file in RT simulations.'})
    use_ramses_acc_rate: bool = field(default=True, 
        metadata={'desc': r' Use RAMSES retion rate field(s) directly (if available) to create the star(s) file in RT simulations.'})
    use_multi_grain: bool = field(default=True, 
        metadata={'desc': r'Do RT with multiple bins if True (if RAMSES output has multiple bins). If False, then dust density is computed using dtogas*rho.'})



    def __repr__(self):  # terminal
        return fancy_repr(self)

    def _repr_html_(self):  # Jupyter
        return _html_repr(self)
    

# ---------------- Structure Params ----------------
@dataclass
class ConfigParams(StrictDataclass):
    ramsesoutput: RamsesOutput= field(default_factory=RamsesOutput)
    radoutput: RadiativeOutput= field(default_factory=RadiativeOutput)
    pymsesrc: Pymsesrc = field(default_factory=Pymsesrc)
    sim: Sim = field(default_factory=Sim)

    def __repr__(self):
        return (
            "PARAMS\n-------\n"
            + fancy_repr(self.inout)
            + "\n\n"
            + fancy_repr(self.pymsesrc)
            + "\n\n"
            + fancy_repr(self.sim)
        )

    def _repr_html_(self):
        body = (
            _html_table(self.inout)
            + "<br/>"
            + _html_table(self.pymsesrc)
            + "<br/>"
            + _html_table(self.sim)
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