import os
import time
import sys
import io
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

logo1_b64 = encode_image_base64(os.path.join(STATIC_DIR, "alma-mater-logo-vector.png"))
logo2_b64 = encode_image_base64(os.path.join(STATIC_DIR, "icelogo.png"))
logo3_b64 = encode_image_base64(os.path.join(STATIC_DIR, "ecogal2.png"))
logo4_b64 = encode_image_base64(os.path.join(STATIC_DIR, "logo-dark.png"))


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
pipe = Pipeline()
cfg = pipe.configparams

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

    status_text = f"RAMSES data will be loaded from:\n{ramses_dir}\n"

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
            <div style="font-size:32px; font-weight:800; display:flex; align-items:center; gap:5px;">
                <span style="display:flex; align-items:center; transform: translateY(7px);">astroMUGS</span>
                <img src="data:image/png;base64,{logo4_b64}" style="height:70px; display:block; transform: translateY(-10px);">
            </div>
            
            <div style="display:flex; gap:20px; align-items:center;">
                <img src="data:image/png;base64,{logo1_b64}" style="height:50px;">
                <img src="data:image/png;base64,{logo2_b64}" style="height:50px;">
                <img src="data:image/png;base64,{logo3_b64}" style="height:50px;">
            </div>
        </div>
        """)


        # # ---------- RAMSES Section ----------
        # gr.Markdown("""
        # # **RAMSES Section**
        # <hr style="border:2px solid #333; margin-top:5px; margin-bottom:20px;">
        # """)

        # ---------- Main Row ----------
        with gr.Row():
            # Main Tabs: RAMSES section + RADMC3D section
            with gr.Tabs() as main_tabs_row:

                # --------- RAMSES Tab ---------
                with gr.Tab("RAMSES Section"):

                    # ---------- Inner Row for RAMSES content ----------
                    with gr.Row():
                        # Left column: input + button + status
                        with gr.Column(scale=1):
                            ramses_dir_input = gr.Textbox(
                                label="Enter the RAMSES output directory (absolute path recommended):",
                                value="ramses_outputs/",
                                lines=1,
                                placeholder="enter path to RAMSES output here/"
                            )

                            start_button = gr.Button("Read RAMSES", variant="primary")

                            status_box = gr.Textbox(
                                label="Status",
                                value="Awaiting input...",
                                lines=4
                            )

                        # Right column: hydro/sink/pymsesrc tabs
                        with gr.Column(scale=2):
                            with gr.Tabs() as ram_tabs:
                                with gr.Tab("hydro_file_descriptor"):
                                    hydro_html = gr.HTML(value="", elem_id="hydro_display")
                                with gr.Tab("sink info"):
                                    sink_html = gr.HTML(value="", elem_id="sink_display")
                                with gr.Tab("pymsesrc"):
                                    pymsesrc_html = gr.HTML(value="", elem_id="pymsesrc_display")

                    # ---------- Hidden states ----------
                    sink_data_state = gr.State()
                    sink_columns_state = gr.State()  # To store column names

                    # ---------- Second row: Plot + Update tabs ----------
                    with gr.Row():
                        with gr.Column(scale=3):
                            with gr.Tabs() as plot_tabs:

                                # Existing Plot Tab
                                with gr.Tab("Plot sinks"):
                                    x_dropdown = gr.Dropdown(choices=[], label="X Column", allow_custom_value=True)
                                    y_dropdown = gr.Dropdown(choices=[], label="Y Column", allow_custom_value=True)
                                    plot_button = gr.Button("Plot", variant="secondary")
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

                                # Update pymsesrc Tab
                                # ---------------- Update pymsesrc Tab ----------------
                                with gr.Tab("Update pymsesrc"):
                                    gr.Markdown("### Pymsesrc Parameters")

                                    py_fields = fields(cfg.pymsesrc)
                                    py_keys = []
                                    py_widgets = []

                                    for f in py_fields:
                                        f_name = f.name
                                        f_default = getattr(cfg.pymsesrc, f_name)
                                        f_desc = f.metadata.get("desc", "")

                                        py_keys.append(f_name)

                                        widget = gr.Checkbox(label=f_desc, value=f_default)
                                        py_widgets.append(widget)

                                    update_pymses_button = gr.Button("Update", variant="primary")
                                    pymses_update_status = gr.Textbox(label="Update Status", value="No updates yet.", lines=2)

                                    def update_pymsesrc_wrapper(*values):
                                        global cfg, pipe

                                        # 1. Update in-memory config
                                        for key, val in zip(py_keys, values):
                                            try:
                                                setattr(cfg.pymsesrc, key, bool(val))
                                            except Exception as e:
                                                return f"Failed to set '{key}': {e}", gr.update(value=None)

                                        # 2. Save pymsesrc
                                        try:
                                            pipe.set_pymsesrc()
                                        except Exception as e:
                                            return f"Updated memory but failed to write file: {e}", gr.update(value=None)

                                        # 3. Reload HTML table
                                        try:
                                            pymses_dict = pipe.read_pymsesrc()
                                            new_html = dict_to_table_html(pymses_dict)
                                            new_html = f"""
                                            <div style="overflow:auto; max-height:400px; border:1px solid #888; border-radius:6px; padding:5px;">
                                                {new_html}
                                            </div>
                                            """
                                        except Exception as e:
                                            return f"Updated but failed to reload table: {e}", gr.update(value=None)

                                        return "✔ Pymsesrc updated successfully.", gr.update(value=new_html)
                                    
                    
                                update_pymses_button.click(
                                    fn=update_pymsesrc_wrapper,
                                    inputs=py_widgets,
                                    outputs=[pymses_update_status, pymsesrc_html]
                                )



                                # Update Parameters Tab
                                with gr.Tab("Update Parameters"):
                                    sim_fields = fields(cfg.sim)
                                    sim_widget_keys = []
                                    sim_widgets = []

                                    gr.Markdown("### Simulation Parameters")

                                    for f in sim_fields:
                                        f_name = f.name
                                        f_type = f.type
                                        f_default = getattr(cfg.sim, f_name)
                                        f_desc = f.metadata.get("desc", "")

                                        sim_widget_keys.append(f_name)

                                        if f_type == bool:
                                            widget = gr.Checkbox(label=f"{f_desc}", value=f_default)
                                        elif f_type in (int, float):
                                            widget = gr.Number(label=f"{f_desc}", value=f_default)
                                        else:
                                            widget = gr.Textbox(label=f"{f_desc}", value=str(f_default))

                                        sim_widgets.append(widget)

                                    update_sim_button = gr.Button("Update", variant="primary")
                                    sim_update_status = gr.Textbox(label="Update Status", value="No updates yet.", lines=2)

                                    def update_sim_params_wrapper(*values):
                                        global cfg, pipe
                                        for key, val in zip(sim_widget_keys, values):
                                            f_decl = next(ff for ff in sim_fields if ff.name == key)
                                            expected = f_decl.type
                                            try:
                                                if expected is bool:
                                                    setattr(cfg.sim, key, bool(val))
                                                elif expected is int:
                                                    setattr(cfg.sim, key, int(val))
                                                elif expected is float:
                                                    setattr(cfg.sim, key, float(val))
                                                else:
                                                    setattr(cfg.sim, key, val)
                                            except Exception as e:
                                                return f"❌ Failed to set '{key}': {e}"
                                        try:
                                            pipe = Pipeline()
                                            pipe.configparams = cfg
                                        except Exception as e:
                                            return f"⚠️ Updated params but pipeline rebuild failed: {e}"
                                        return "✔ Simulation parameters updated successfully."

                                    update_sim_button.click(
                                        fn=update_sim_params_wrapper,
                                        inputs=sim_widgets,
                                        outputs=sim_update_status
                                    )


                # --------- RADMC3D Tab (empty) ---------
                with gr.Tab("RADMC3D Section"):
                    # ---------- Inner Row for RAMSES content ----------
                    with gr.Row():
                        # Left column: input + button + status
                        with gr.Column(scale=1):
                            radmc3d_dir_input = gr.Textbox(
                                label="Enter the RADMC3D output directory (absolute path recommended):",
                                value="radmc3d_outputs/",
                                lines=1,
                                placeholder="enter path to RADMC3D output here/"
                            )

                            set_radmc3d_button = gr.Button("Set RADMC3D directory", variant="primary")


                            radmc_status_box = gr.Textbox(
                                label="Status",
                                value="Awaiting input...",
                                lines=4
                            )

                # --------- POLARIS Tab (empty) ---------
                with gr.Tab("POLARIS Section"):
                    gr.Markdown("*(This tab is currently empty and reserved for POLARIS controls.)*")

                # --------- PIPLINE Tab ---------
                with gr.Tab("PIPELINE"):
                    gr.Markdown("### Run RAMSES → RADMC3D Conversion")
                    convert_button = gr.Button("Convert to RADMC3D", variant="primary")
                    convert_log = gr.Textbox(
                        label="Live Conversion Log",
                        value="Logs will appear here...",
                        lines=20,
                        interactive=False
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

        # ---------- Callback for Read RAMSES button ----------
        def on_set_radmc3d(radmc_dir):
            try:
                info = pipe.set_radmc3d_dir(radmc_dir)
                return info
            except Exception as e:
                return f"Error: {e}"
            
        set_radmc3d_button.click(
            fn=on_set_radmc3d,
            inputs=[radmc3d_dir_input],
            outputs=[radmc_status_box]
        )

        # ---------- Convert to RADMC3D callback ----------
        def on_convert():
            try:
                pipe.convert_rasmes2radmc()
                return "✔ Conversion complete."
            except Exception as e:
                return f"❌ Conversion failed: {e}"

        def on_convert_stream():
            """
            Streams output line-by-line to the Gradio UI during conversion.
            """

            # Capture stdout
            buffer = io.StringIO()
            old_stdout = sys.stdout
            sys.stdout = buffer

            yield "Starting RAMSES → RADMC3D conversion...\n"

            try:
                # Surround the conversion with our own messaging
                print("Initializing conversion...")
                yield buffer.getvalue(); buffer.truncate(0); buffer.seek(0)

                pipe.convert_rasmes2radmc()

                # Capture the final print statements
                yield buffer.getvalue(); buffer.truncate(0); buffer.seek(0)

                yield "✔ Conversion complete.\n"

            except Exception as e:
                yield buffer.getvalue()

                yield f"❌ Conversion failed: {e}\n"

            finally:
                # Restore stdout no matter what
                sys.stdout = old_stdout

        convert_button.click(
            fn=on_convert_stream,
            inputs=[],
            outputs=convert_log
        )

        # ---------- Update dropdowns when column names are ready ----------
        def update_dropdowns(columns):
            if not columns:
                return gr.update(choices=[], value=None), gr.update(choices=[], value=None)
            return (
                gr.update(choices=columns, value=columns[0]),
                gr.update(choices=columns, value=columns[1] if len(columns) > 1 else columns[0])
            )

        sink_columns_state.change(
            fn=update_dropdowns,
            inputs=[sink_columns_state],
            outputs=[x_dropdown, y_dropdown]
        )

    demo.launch(theme=gr.themes.Citrus())
    return demo



# -----------------------------
# MAIN
# -----------------------------
if __name__ == "__main__":
    demo = launch_interface()
    demo.launch(share=True)
