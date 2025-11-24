import os
import inspect
from dataclasses import fields

import base64
import gradio as gr
from gradio import Dropdown

from radprocess.utils.config import ConfigParams
from radprocess.pipeline.Pipeline import Pipeline
from radprocess.plotting.plot import plot_sink_columns

# -----------------------------
# Paths
# -----------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
STATIC_DIR = os.path.join(BASE_DIR, "_static")

def encode_image_base64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")

logo1_b64 = encode_image_base64(os.path.join(STATIC_DIR, "icelogo.png"))
logo2_b64 = encode_image_base64(os.path.join(STATIC_DIR, "ecogal2.png"))


def dict_to_table_html(d, level=0):
    """Recursively converts a dict into a nested HTML table with colored keys and values."""
    html_str = '<table style="border-collapse: collapse; font-family: monospace; width:100%;">'
    
    for key, val in d.items():
        html_str += "<tr>"
        # Key column (blue, indented by level)
        html_str += f'<td style="padding:4px; font-weight:bold; color:#2a52be; vertical-align:top;">{"&nbsp;"*4*level}{key}</td>'
        
        # Value column
        if isinstance(val, dict):
            html_str += f'<td style="padding:4px;">{dict_to_table_html(val, level+1)}</td>'
        elif isinstance(val, list):
            html_str += '<td style="padding:4px;">'
            for item in val:
                if isinstance(item, dict):
                    html_str += dict_to_table_html(item, level+1)
                else:
                    html_str += f'{item}<br>'
            html_str += '</td>'
        else:
            html_str += f'<td style="padding:4px; color:#008000;">{val}</td>'
        
        html_str += "</tr>"
    
    html_str += "</table>"
    return html_str

def sink_text_to_table(text):
    lines = text.strip().split("\n")
    if not lines:
        return "<div>No sink info available</div>"

    # First line is header (column names)
    header = lines[0].split()
    data_rows = [line.split() for line in lines[1:]]

    table_html = f"""
    <div style="overflow:auto; max-height:400px; border:1px solid #828181; border-radius:6px; padding:5px;">
        <table style="border-collapse: collapse; width:100%; font-family: monospace;">
            <thead>
                <tr style="background-color:#737373; color:white; position: sticky; top: 0;">
                    {''.join(f'<th style="padding:6px; border:1px solid #F2F2F2;">{col}</th>' for col in header)}
                </tr>
            </thead>
            <tbody>
                {''.join(
                    f'<tr style="background-color:{"#2E2E2E" if i%2==0 else "black"};">'
                    + ''.join(f'<td style="padding:4px; border:1px solid #ccc;">{cell}</td>' for cell in row)
                    + '</tr>'
                    for i, row in enumerate(data_rows)
                )}
            </tbody>
        </table>
    </div>
    """
    return table_html


# -----------------------------
# Global config and pipeline
# -----------------------------
cfg = ConfigParams()
pipe = None

# -----------------------------
# Core update function
# -----------------------------
def start_pipeline_display(ramses_dir):
    global cfg, pipe
    cfg.ramsesoutput.ramses_output_dir = ramses_dir
    pipe = Pipeline()
    pipe.configparams = cfg

    hydro_text = pipe.read_hydro_descriptor()
    sink_text = pipe.read_sink_info()
    pymses_dict = pipe.read_pymsesrc()

    # Convert sink text to rows for plotting
    lines = sink_text.strip().split("\n")
    header = lines[0].split() if lines else []
    sink_rows = [dict(zip(header, line.split())) for line in lines[1:]] if len(lines) > 1 else []

    status_text = f"RAMSES directory set to:\n{ramses_dir}\n"

    hydro_html = f"""
    <div style="overflow:auto; max-height:400px; border:1px solid #888; border-radius:6px; padding:5px;">
        <pre style="font-family: monospace; margin:0;">{hydro_text}</pre>
    </div>
    """

    sink_html = sink_text_to_table(sink_text)

    pymsesrc_html = f"""
    <div style="overflow:auto; max-height:400px; border:1px solid #888; border-radius:6px; padding:5px;">
        {dict_to_table_html(pymses_dict)}
    </div>
    """

    # Return sink_rows as the 5th element
    return status_text, hydro_html, sink_html, pymsesrc_html, sink_rows




# -----------------------------
# Gradio Layout
# -----------------------------
def launch_interface():
    with gr.Blocks(title="radprocess Interface") as demo:

        # ---------- Banner ----------
        gr.HTML(f"""
        <div style="
            display:flex; 
            justify-content:space-between; 
            align-items:center;
            padding: 15px 20px;
            background: #1a1a1a;
            color:white;
            font-family: sans-serif;
            border-bottom: 1px solid #666;
            position:sticky;
            top:0;
            z-index:100;
        ">
            <div style="font-size:32px; font-weight:800;">radprocess</div>
            <div style="display:flex; gap:20px;">
                <img src="data:image/png;base64,{logo1_b64}" style="height:50px;">
                <img src="data:image/png;base64,{logo2_b64}" style="height:50px;">
            </div>
        </div>
        """)

        # ---------- RAMSES Section ----------
        gr.Markdown("""
        # **RAMSES Section**
        <hr style="border:2px solid #333; margin-top:5px; margin-bottom:20px;">
        """)

        # ---------- Main Row ----------
        with gr.Row():
            # Left column: input + button + status
            with gr.Column(scale=1):
                ramses_dir_input = gr.Textbox(
                    label="Enter the RAMSES output directory (absolute path recommended):",
                    value="ramses_outputs/",
                    lines=1,
                    placeholder="/absolute/path/to/ramses/output/"
                )

                start_button = gr.Button("Read RAMSES", variant="primary")

                status_box = gr.Textbox(
                    label="Status",
                    value="Awaiting input...",
                    lines=4
                )

            # Right column: tabs for hydro/sink/pymsesrc
            with gr.Column(scale=2):
                with gr.Tabs() as main_tabs:
                    with gr.Tab("hydro_file_descriptor"):
                        hydro_html = gr.HTML(value="", elem_id="hydro_display")
                    with gr.Tab("sink info"):
                        sink_html = gr.HTML(value="", elem_id="sink_display")
                    with gr.Tab("pymsesrc"):
                        pymsesrc_html = gr.HTML(value="", elem_id="pymsesrc_display")

        # ---------- Hidden states ----------
        sink_data_state = gr.State()
        sink_columns_state = gr.State()  # To store column names

        # ---------- Second row ----------
        with gr.Row():
            with gr.Column(scale=3):
                with gr.Tabs() as plot_tabs:

                    # ----------------- Existing Plot Tab -----------------
                    with gr.Tab("Plot sinks"):
                        # Dropdowns first
                        x_dropdown = gr.Dropdown(choices=[], label="X Column", allow_custom_value=True)
                        y_dropdown = gr.Dropdown(choices=[], label="Y Column", allow_custom_value=True)

                        # Plot button below dropdowns
                        plot_button = gr.Button("Plot", variant="secondary")

                        # Plot output
                        plot_output = gr.Plot()

                        def generate_sink_plot(sink_rows, x_col, y_col):
                            if not sink_rows or not x_col or not y_col:
                                return None
                            fig = plot_sink_columns(sink_rows, x_col, y_col)
                            return fig

                        plot_button.click(
                            fn=generate_sink_plot,
                            inputs=[sink_data_state, x_dropdown, y_dropdown],
                            outputs=[plot_output]
                        )

                    # ----------------- NEW EMPTY TABS -----------------
                    with gr.Tab("Update pymsesrc"):
                        gr.Markdown("*(This tab will be used to configure parameter set A.)*")

                    with gr.Tab("Update Parameters"):

                        # --- Build dynamic UI based on Sim dataclass ---
                        sim_fields = fields(cfg.sim)

                        # Use an ordered list to keep stable mapping between fields and values
                        sim_widget_keys = []
                        sim_widgets = []

                        gr.Markdown("### Simulation Parameters")

                        for f in sim_fields:
                            f_name = f.name
                            f_type = f.type
                            f_default = getattr(cfg.sim, f_name)
                            f_desc = f.metadata.get("desc", "")

                            sim_widget_keys.append(f_name)

                            # Boolean -> Checkbox (nice native widget)
                            if f_type == bool:
                                widget = gr.Checkbox(label=f"{f_desc}", value=f_default)

                            # Numeric -> Number
                            elif f_type in (int, float):
                                # gr.Number supports both int/float; keep value as f_default
                                widget = gr.Number(label=f"{f_desc}", value=f_default)

                            else:
                                # Fallback: text field (string/untyped)
                                widget = gr.Textbox(label=f"{f_desc}", value=str(f_default))

                            sim_widgets.append(widget)

                        # Status output and Update button
                        update_sim_button = gr.Button("Update", variant="primary")
                        sim_update_status = gr.Textbox(label="Update Status", value="No updates yet.", lines=2)

                        # Wrapper that receives positional args in the same order as sim_widgets list
                        def update_sim_params_wrapper(*values):
                            """
                            values is a tuple with one element per widget in sim_widgets,
                            in the same order as sim_widget_keys.
                            """
                            global cfg, pipe
                            # Map values to field names and set attributes (StrictDataclass will validate)
                            for key, val in zip(sim_widget_keys, values):
                                # Convert Textbox fallback strings to original type if needed
                                f_decl = next(ff for ff in sim_fields if ff.name == key)
                                expected = f_decl.type
                                try:
                                    if expected is bool:
                                        # Checkbox already returns bool
                                        setattr(cfg.sim, key, bool(val))
                                    elif expected is int:
                                        setattr(cfg.sim, key, int(val))
                                    elif expected is float:
                                        setattr(cfg.sim, key, float(val))
                                    else:
                                        # fallback: store string
                                        setattr(cfg.sim, key, val)
                                except Exception as e:
                                    return f"❌ Failed to set '{key}': {e}"

                            # Rebuild pipeline to apply new config (optional)
                            try:
                                pipe = Pipeline()
                                pipe.configparams = cfg
                            except Exception as e:
                                return f"⚠️ Updated params but pipeline rebuild failed: {e}"

                            return "✔ Simulation parameters updated successfully."

                        # Connect button: pass the list of components (NOT a dict)
                        update_sim_button.click(
                            fn=update_sim_params_wrapper,
                            inputs=sim_widgets,   # list of components
                            outputs=sim_update_status
                        )




        # ---------- Callback for Read RAMSES button ----------
        def on_read_ramses(ramses_dir):
            status, hydro_html_val, sink_html_val, pymsesrc_html_val, sink_rows = start_pipeline_display(ramses_dir)

            # Extract columns from the first row of sink_rows
            if sink_rows:
                columns = list(sink_rows[0].keys())
            else:
                columns = []

            return (
                status,
                hydro_html_val,
                sink_html_val,
                pymsesrc_html_val,
                sink_rows,
                columns  # Update sink_columns_state
            )

        start_button.click(
            fn=on_read_ramses,
            inputs=[ramses_dir_input],
            outputs=[
                status_box,
                hydro_html,
                sink_html,
                pymsesrc_html,
                sink_data_state,
                sink_columns_state
            ]
        )

        # ---------- Update dropdowns when column names are ready ----------
        def update_dropdowns(columns):
            if not columns:
                return gr.Dropdown.update(choices=[], value=None), gr.Dropdown.update(choices=[], value=None)
            return (
                gr.Dropdown.update(choices=columns, value=columns[0]),
                gr.Dropdown.update(choices=columns, value=columns[1] if len(columns) > 1 else columns[0])
            )

        sink_columns_state.change(
            fn=update_dropdowns,
            inputs=[sink_columns_state],
            outputs=[x_dropdown, y_dropdown]
        )

    return demo



# -----------------------------
# MAIN
# -----------------------------
if __name__ == "__main__":
    demo = launch_interface()
    demo.launch()
