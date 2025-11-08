from dataclasses import dataclass, fields, is_dataclass, field
import html
from typing import Literal, Any

import numpy as np


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

def html_table_from_columns(columns, rows, title="Table"):
    """Create a collapsible HTML table from arbitrary columns + rows."""
    header_cells = "".join(
        f'<th style="padding-right:20px; text-align:left;">{col}</th>'
        for col in columns
    )

    body_rows = ""
    for row in rows:
        if isinstance(row, dict):
            values = [row[c] for c in columns]
        else:
            values = row

        cells = "".join(
            f'<td style="text-align:left;">{html_value(v)}</td>' for v in values
        )
        body_rows += f"<tr>{cells}</tr>\n"

    return f"""
    <details style="margin:8px 0; padding:6px; border:1px solid #ccc; border-radius:6px;">
      <summary style="font-size:16px; font-weight:bold; cursor:pointer;">{title}</summary>
      <div style="padding:10px;">

        <table style="border-collapse: collapse; text-align:left;">
            <thead>
                <tr style="border-bottom:1px solid #888;">
                    {header_cells}
                </tr>
            </thead>
            <tbody>
                {body_rows}
            </tbody>
        </table>

      </div>
    </details>

    <script>
    if (window.MathJax) {{
        MathJax.typesetPromise();
    }}
    </script>
    """


# =======================================================
# PARAMETER DATACLASSES
# =======================================================
@dataclass
class SinkInfo:
    columns: list
    rows: list          # string-based rows (from .info)
    num_sinks: int
    data: np.ndarray    # numeric data from CSV

    def _repr_html_(self):
        title = f"Sink File ({self.num_sinks} sink{'s' if self.num_sinks != 1 else ''})"
        return html_table_from_columns(self.columns, self.rows, title=title)

    def __repr__(self):
        return f"SinkTable(num_sinks={self.num_sinks}, columns={self.columns})"