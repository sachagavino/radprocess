import os
import time
import sys
import io
from dataclasses import fields
import threading
import time

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
def start_pipeline_display(ramses_dir, cfg, pipe):
    # Update config
    cfg.ramsesoutput.ramses_output_dir = ramses_dir
    # Rebuild pipeline
    pipe = Pipeline()
    pipe.configparams = cfg

    status_text = f"RAMSES data will be loaded from:\n{ramses_dir}\n"

    # AUTO-WRITE ~/.pymses/pymsesrc
    try:
        pipe.set_pymsesrc()
    except Exception as e:
        status_text += f"\n⚠ pymsesrc write failed: {e}"

    hydro_text = pipe.read_hydro_descriptor()
    other_hydro_text = pipe.read_other_file_descriptor()
    sink_text = pipe.read_sink_info()
    pymses_dict = pipe.read_pymsesrc()

    lines = sink_text.strip().split("\n")
    header = lines[0].split() if lines else []
    sink_rows = [dict(zip(header, line.split())) for line in lines[1:]] if len(lines) > 1 else []

    hydro_html = f"""
    <div style="overflow:auto; max-height:400px; border:1px solid #888; border-radius:6px; padding:5px;">
        <pre style="font-family: monospace; margin:0;">{hydro_text}</pre>
    </div>
    """

    other_hydro_html = f"""
    <div style="overflow:auto; max-height:400px; border:1px solid #888; border-radius:6px; padding:5px;">
        <pre style="font-family: monospace; margin:0;">{other_hydro_text}</pre>
    </div>
    """

    sink_html = sink_text_to_table(sink_text)
    pymsesrc_html = f"""
    <div style="overflow:auto; max-height:400px; border:1px solid #888; border-radius:6px; padding:5px;">
        {dict_to_table_html(pymses_dict)}
    </div>
    """

    return status_text, hydro_html, other_hydro_html, sink_html, pymsesrc_html, sink_rows, cfg, pipe


# -----------------------------
# Gradio Layout
# -----------------------------
def launch_interface():
    # 1. Define all components first
    status_box = gr.Textbox(...)
    hydro_html = gr.HTML(value="")
    other_hydro_html = gr.HTML(value="")
    sink_html = gr.HTML(value="")
    pymsesrc_html = gr.HTML(value="")



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

        # ---------- Convert callback ----------
        def on_convert():
            try:
                pipe.load_ramses()
                #return "✔ Conversion complete."
            except Exception as e:
                return f"❌ Conversion failed: {e}"


        ##-- STREAM LOADING RAMSES --
        def on_convert_stream(cfg, pipe):
            buffer = io.StringIO()
            old_stdout = sys.stdout
            sys.stdout = buffer

            yield "Starting RAMSES → AMR grid conversion...\n"

            try:
                print("Initializing conversion...")
                yield buffer.getvalue(); buffer.truncate(0); buffer.seek(0)

                # IMPORTANT: use updated pipeline
                pipe.configparams = cfg
                pipe.load_ramses()

                yield buffer.getvalue(); buffer.truncate(0); buffer.seek(0)

            except Exception as e:
                yield buffer.getvalue()
                yield f"❌ Conversion failed: {e}\n"

            finally:
                sys.stdout = old_stdout

        ##-- STREAM CONVERTINT RAMSES TO RADMC3D --
        def on_convert_radmc_stream(cfg, pipe):
            buffer = io.StringIO()
            old_stdout = sys.stdout
            sys.stdout = buffer

            yield "Starting AMR → RADMC3D conversion...\n"

            try:
                pipe.configparams = cfg

                print("Initializing RADMC3D conversion...")
                yield buffer.getvalue(); buffer.truncate(0); buffer.seek(0)

                pipe.convert_to_radmc()

                yield buffer.getvalue(); buffer.truncate(0); buffer.seek(0)

            except Exception as e:
                yield buffer.getvalue()
                yield f"❌ RADMC3D conversion failed: {e}\n"

            finally:
                sys.stdout = old_stdout



        # ---------- Main Row ----------
        with gr.Row():
            # Main Tabs: RAMSES section
            with gr.Tabs() as main_tabs_row:

                # --------- RAMSES Tab ---------
                with gr.Tab("RAMSES Section"):

                    # ---------- Hidden states ----------
                    sink_data_state = gr.State([])
                    sink_columns_state = gr.State([])
                    cfg_state = gr.State(cfg)
                    pipe_state = gr.State(pipe)

                    # ---------- Inner Row for RAMSES content ----------
                    with gr.Row():

                        # ===== LEFT COLUMN =====
                        with gr.Column(scale=1):

                            ramses_dir_input = gr.Textbox(
                                label="Enter the RAMSES output directory (absolute path recommended):",
                                value=None,
                                lines=1,
                                placeholder="enter path to RAMSES output here/"
                            )

                            start_button = gr.Button("Load RAMSES output", variant="primary")

                            status_box = gr.Textbox(
                                label="Status",
                                value="Awaiting input...",
                                lines=4
                            )

                        # ===== RIGHT COLUMN =====
                        with gr.Column(scale=2):

                            with gr.Tabs() as ram_tabs:

                                # --- hydro descriptor ---
                                with gr.Tab("hydro_file_descriptor"):
                                    hydro_html = gr.HTML(value="", elem_id="hydro_display")

                                # --- hydro descriptor ---
                                with gr.Tab("other_file_descriptor"):
                                    other_hydro_html = gr.HTML(value="", elem_id="other_hydro_display")


                                # --- pymsesrc ---
                                with gr.Tab("pymsesrc"):
                                    pymsesrc_html = gr.HTML(value="", elem_id="pymsesrc_display")

                                # --- sink info ---
                                with gr.Tab("sink info"):
                                    sink_html = gr.HTML(value="", elem_id="sink_display")

                                # --- NEW: Plot sinks ---
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





                    # ---------- Second row: Plot + Update tabs ----------
                    with gr.Row():
                        with gr.Column(scale=3):
                            with gr.Tabs() as plot_tabs:

                                # # Existing Plot Tab
                                # with gr.Tab("Plot sinks"):
                                #     x_dropdown = gr.Dropdown(choices=[], label="X Column", allow_custom_value=True)
                                #     y_dropdown = gr.Dropdown(choices=[], label="Y Column", allow_custom_value=True)
                                #     plot_button = gr.Button("Plot", variant="secondary")
                                #     plot_output = gr.Plot()

                                #     def generate_sink_plot(sink_rows, x_col, y_col):
                                #         if not sink_rows or not x_col or not y_col:
                                #             return None
                                #         fig = plot_sink_columns(sink_rows, x_col, y_col)
                                #         return fig

                                #     plot_button.click(
                                #         fn=generate_sink_plot,
                                #         inputs=[sink_data_state, x_dropdown, y_dropdown],
                                #         outputs=[plot_output]
                                #     )

                                # Update amr source Tab
                                # ---------------- Update amr source Tab ----------------
                                with gr.Tab("Set amr source"):
                                    gr.Markdown("### Select here which variables you want to include in the conversion")

                                    py_fields = fields(cfg.amrsource)
                                    py_keys = []
                                    py_widgets = []

                                    for f in py_fields:
                                        f_name = f.name
                                        f_default = getattr(cfg.amrsource, f_name)
                                        f_desc = f.metadata.get("desc", "")

                                        py_keys.append(f_name)

                                        widget = gr.Checkbox(label=f_desc, value=f_default)
                                        py_widgets.append(widget)

                                    update_amrsource_button = gr.Button("Update", variant="primary")
                                    amrsource_update_status = gr.Textbox(label="Update Status", value="No updates yet.", lines=2)


                                    def update_amrsource_wrapper(*values_and_states):
                                        *values, cfg, pipe = values_and_states
                                        # 1. Update in-memory config
                                        for key, val in zip(py_keys, values):
                                            try:
                                                setattr(cfg.amrsource, key, bool(val))
                                            except Exception as e:
                                                return f"❌ Failed to set '{key}': {e}"

                                        try:
                                            pipe = Pipeline()
                                            pipe.configparams = cfg
                                        except Exception as e:
                                            return f"⚠️ Updated amr source but pipeline rebuild failed: {e}"
                                        return "✔ AMR source updated successfully.", cfg, pipe

                                update_amrsource_button.click(
                                    fn=update_amrsource_wrapper,
                                    inputs=py_widgets + [cfg_state, pipe_state],
                                    outputs=[amrsource_update_status, cfg_state, pipe_state]
                                )

                                # Update Parameters Tab
                                with gr.Tab("Set Parameters"):
                                    gr.Markdown("### Simulation Parameters")

                                    sim_fields = fields(cfg.sim)
                                    sim_widget_keys = []
                                    sim_widgets = []

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

                                    def update_sim_params_wrapper(*values_and_states):
                                        # values_and_states = [val1, val2, ..., cfg, pipe]
                                        *values, cfg, pipe = values_and_states

                                        # 1) Update in-memory config object
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

                                        # 2) Recreate pipeline object with updated config
                                        try:
                                            pipe = Pipeline()
                                            pipe.configparams = cfg
                                        except Exception as e:
                                            return f"⚠️ Updated params but pipeline rebuild failed: {e}"

                                        # 3) Return status + updated states so Gradio stores them
                                        return "✔ Simulation parameters updated successfully.", cfg, pipe


                                    update_sim_button.click(
                                        fn=update_sim_params_wrapper,
                                        inputs=sim_widgets + [cfg_state, pipe_state],
                                        outputs=[sim_update_status, cfg_state, pipe_state]
                                    )




                # --------- RADMC3D Tab (empty) ---------
                with gr.Tab("RADMC3D Section"):
                    with gr.Tabs() as radmc_tabs_row:
                        # --------- Working directory Tab ---------
                        with gr.Tab("radmc3d.inp"):
                            print('nothing')
                        with gr.Tab("wavelength_microns.inp"):
                            print('nothing')
                        with gr.Tab("stars.inp"):
                            print('nothing')

                # --------- POLARIS Tab (empty) ---------
                with gr.Tab("POLARIS Section"):
                    # ---------- Inner Row for POLARIS content ----------
                    with gr.Row():
                        # Left column: input + button + status
                        with gr.Column(scale=1):
                            polaris_dir_input = gr.Textbox(
                                label="Enter the POLARIS output directory (absolute path recommended):",
                                value="polaris_outputs/",
                                lines=1,
                                placeholder="enter path to POLARIS output here/"
                            )

                            # set_polaris_button = gr.Button("Set POLARIS directory", variant="primary")


                            # polaris_status_box = gr.Textbox(
                            #     label="Status",
                            #     value="Awaiting input...",
                            #     lines=4
                            # )

                # --------- PIPELINE Tab ---------
                with gr.Tab("PIPELINE Section"):
                    with gr.Tabs() as pipeline_tabs_row:
                        # --------- Working directory Tab ---------
                        with gr.Tab("Working directory"):

                            # --------- Working directory (top of PIPELINE tab) ---------
                            with gr.Row():
                                with gr.Column(scale=1):

                                    workdir_input = gr.Textbox(
                                        label="Enter the working directory (absolute path recommanded):",
                                        value=None,
                                        lines=1,
                                        placeholder="e.g. model_output/"
                                    )

                                    set_workdir_button = gr.Button("Set working directory", variant="primary")

                                    workdir_status_box = gr.Textbox(
                                        label="Working directory status",
                                        value="Awaiting input...",
                                        lines=3
                                    )




                        # --------- MAIN Pipeline Tab ---------
                        with gr.Tab("PIPELINE"): 
                            with gr.Tab("Interactive"):  
                                #DEFINE THE ACTIONS HERE.
                                BUTTON_VISIBILITY = {
                                    "Create RAMSES grid": True,
                                    "Convert to RADMC3D": False,
                                    "Convert to POLARIS": False,
                                    "merge opacities (POLARIS to RADMC3D)": False,
                                    "Run MCTHERM (RADMC3D)": False,
                                    "Merge temperatures": False,
                                    "Make image": False,
                                }

                                pipeline_actions = list(BUTTON_VISIBILITY.keys())
                                MAX_STEPS = 10

                                step_count = gr.State(0)                           

                                with gr.Row():

                                    # ===== LEFT COLUMN =====
                                    with gr.Column(scale=1.5):

                                        gr.Markdown("### Build the pipeline by selecting actions below.")

                                        rows = []
                                        dropdowns = []
                                        buttons = []

                                        for i in range(MAX_STEPS):

                                            # 🔹 Group keeps dropdown & button together
                                            with gr.Group(visible=False) as row:

                                                # --- OPEN BOX ---
                                                gr.HTML(
                                                    "<div style='border:1px solid #555;"
                                                    "padding:10px;"
                                                    "border-radius:8px;"
                                                    "margin-bottom:8px;'>"
                                                )

                                                with gr.Row():

                                                    dd = gr.Dropdown(
                                                        choices=pipeline_actions,
                                                        label=f"Action {i+1}",
                                                        value=None,
                                                        interactive=True
                                                    )

                                                    btn = gr.Button(
                                                        "Create RAMSES grid",
                                                        visible=False
                                                    )

                                                # --- CLOSE BOX ---
                                                gr.HTML("</div>")

                                            rows.append(row)
                                            dropdowns.append(dd)
                                            buttons.append(btn)

                                        add_step_button = gr.Button("➕ Add action", variant="huggingface")

                                    # ===== RIGHT COLUMN =====
                                    with gr.Column(scale=1.5):
                                        convert_log = gr.Textbox(
                                            label="Live Conversion Log",
                                            value="Logs will appear here...",
                                            lines=20,
                                            interactive=False
                                        )


                                # === increment step count ===
                                def add_step(n):
                                    if n < MAX_STEPS:
                                        n += 1
                                    return n

                                add_step_button.click(
                                    fn=add_step,
                                    inputs=[step_count],
                                    outputs=[step_count]
                                )


                                # === show rows up to step_count ===
                                def update_visibility(n):
                                    return [gr.update(visible=(i < n)) for i in range(MAX_STEPS)]

                                step_count.change(
                                    fn=update_visibility,
                                    inputs=[step_count],
                                    outputs=rows
                                )


                                # === toggle button per dropdown (lambda binding FIXED) ===
                                def update_action_button(action):
                                    if action == "Create RAMSES grid":
                                        return gr.update(visible=True, value="Create RAMSES grid")
                                    elif action == "Convert to RADMC3D":
                                        return gr.update(visible=True, value="Convert to RADMC3D")
                                    else:
                                        return gr.update(visible=False)

                            for dd, btn in zip(dropdowns, buttons):

                                dd.change(
                                    fn=update_action_button,
                                    inputs=[dd],
                                    outputs=[btn],
                                )

                                def make_dispatch():

                                    def dispatch(action, cfg, pipe):
                                        buffer = io.StringIO()
                                        old_stdout = sys.stdout
                                        sys.stdout = buffer

                                        done = False
                                        error = None

                                        full_log = ""

                                        def run_job():
                                            nonlocal done, error
                                            try:
                                                pipe.configparams = cfg

                                                if action == "Create RAMSES grid":
                                                    pipe.load_ramses()

                                                elif action == "Convert to RADMC3D":
                                                    pipe.convert_to_radmc()

                                                else:
                                                    raise RuntimeError("Unknown action.")

                                            except Exception as e:
                                                error = e
                                            finally:
                                                done = True

                                        try:
                                            if action == "Create RAMSES grid":
                                                full_log += "Starting RAMSES → AMR grid conversion...\n"
                                                yield full_log
                                            elif action == "Convert to RADMC3D":
                                                full_log += "Starting AMR → RADMC3D conversion...\n"
                                                yield full_log

                                            t = threading.Thread(target=run_job)
                                            t.start()

                                            while not done:
                                                time.sleep(0.2)
                                                txt = buffer.getvalue()
                                                if txt:
                                                    full_log += txt
                                                    yield full_log
                                                    buffer.truncate(0)
                                                    buffer.seek(0)

                                            txt = buffer.getvalue()
                                            if txt:
                                                full_log += txt
                                                yield full_log

                                            if error is not None:
                                                full_log += f"❌ Failed: {error}\n"
                                                yield full_log
                                            else:
                                                if action == "Create RAMSES grid":
                                                    full_log += "✔ RAMSES grid created successfully.\n"
                                                    yield full_log
                                                elif action == "Convert to RADMC3D":
                                                    full_log += "✔ RADMC3D conversion finished.\n"
                                                    yield full_log

                                        finally:
                                            sys.stdout = old_stdout

                                    return dispatch

                                # THIS is an important part
                                btn.click(
                                    fn=make_dispatch(),
                                    inputs=[dd, cfg_state, pipe_state],
                                    outputs=convert_log,
                                )





                            with gr.Tab("Jobs"):
                                print('empty:')

                        # --------- MAIN Results Tab ---------
                        with gr.Tab("Pipeline Results"):
                            print('nothing here') 



        def on_read_ramses(ramses_dir, cfg, pipe):
            status_text, hydro_html_val, other_hydro_html_val, sink_html_val, pymsesrc_html_val, sink_rows, cfg, pipe = start_pipeline_display(ramses_dir, cfg, pipe)
            # Extract columns for dropdowns
            columns = list(sink_rows[0].keys()) if sink_rows else []
            return status_text, hydro_html_val, other_hydro_html_val, sink_html_val, pymsesrc_html_val, sink_rows, columns, cfg, pipe

        # After all components are defined
        start_button.click(
            fn=on_read_ramses,
            inputs=[ramses_dir_input, cfg_state, pipe_state],
            outputs=[
                status_box,
                hydro_html,
                other_hydro_html,
                sink_html,
                pymsesrc_html,
                sink_data_state,
                sink_columns_state,
                cfg_state,
                pipe_state
            ]
        )





        # # ---------- Callback for Read RADMC3D button ----------
        # def on_set_radmc3d(radmc_dir, cfg, pipe):
        #     try:
        #         # ensure pipeline uses the provided cfg
        #         pipe.configparams = cfg
        #         info = pipe.set_radmc3d_dir(radmc_dir)
        #         return info, cfg, pipe
        #     except Exception as e:
        #         return f"Error: {e}", cfg, pipe

        # set_radmc3d_button.click(
        #     fn=on_set_radmc3d,
        #     inputs=[radmc3d_dir_input, cfg_state, pipe_state],
        #     outputs=[radmc_status_box, cfg_state, pipe_state]
        # )

        # ---------- Callback for Pipeline output ----------
        def on_set_workdir(workdir, cfg, pipe):
            try:
                pipe.configparams = cfg
                msg = pipe.set_working_dir(workdir)
                return msg, cfg, pipe
            except Exception as e:
                return f"Error: {e}", cfg, pipe

        set_workdir_button.click(
            fn=on_set_workdir,
            inputs=[workdir_input, cfg_state, pipe_state],
            outputs=[workdir_status_box, cfg_state, pipe_state]
        )




        # # ---------- Callback for Read POLARIS button ----------
        # def on_set_polaris(polaris_dir, cfg, pipe):
        #     try:
        #         # ensure pipeline uses the provided cfg
        #         pipe.configparams = cfg
        #         info = pipe.set_polaris_dir(polaris_dir)
        #         return info, cfg, pipe
        #     except Exception as e:
        #         return f"Error: {e}", cfg, pipe

        # set_polaris_button.click(
        #     fn=on_set_polaris,
        #     inputs=[polaris_dir_input, cfg_state, pipe_state],
        #     outputs=[polaris_status_box, cfg_state, pipe_state]
        # )



        # create_amr_button.click(
        #     fn=on_convert_stream,
        #     inputs=[cfg_state, pipe_state],
        #     outputs=convert_log
        # )


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

    demo.launch(theme=gr.themes.Citrus(), share=True)
    return demo



# -----------------------------
# MAIN
# -----------------------------
if __name__ == "__main__":
    demo = launch_interface()
    demo.launch()




        # # ---------- Callback for Read RAMSES button ----------
        # def on_read_ramses(ramses_dir, cfg, pipe):
        #     # Update config
        #     cfg.ramsesoutput.ramses_output_dir = ramses_dir
        #     # Rebuild pipeline
        #     pipe = Pipeline()
        #     pipe.configparams = cfg

        #     status_text = f"RAMSES data will be loaded from:\n{ramses_dir}\n"

        #     # AUTO-WRITE ~/.pymses/pymsesrc
        #     try:
        #         pipe.set_pymsesrc()
        #     except Exception as e:
        #         status_text += f"\n⚠ pymsesrc write failed: {e}"

        #     hydro_text = pipe.read_hydro_descriptor()
        #     sink_text = pipe.read_sink_info()
        #     pymses_dict = pipe.read_pymsesrc()

        #     # Convert sink text to rows for plotting
        #     lines = sink_text.strip().split("\n")
        #     header = lines[0].split() if lines else []
        #     sink_rows = [dict(zip(header, line.split())) for line in lines[1:]] if len(lines) > 1 else []

        #     hydro_html = f"""
        #     <div style="overflow:auto; max-height:400px; border:1px solid #888; border-radius:6px; padding:5px;">
        #         <pre style="font-family: monospace; margin:0;">{hydro_text}</pre>
        #     </div>
        #     """
        #     sink_html = sink_text_to_table(sink_text)
        #     pymsesrc_html = f"""
        #     <div style="overflow:auto; max-height:400px; border:1px solid #888; border-radius:6px; padding:5px;">
        #         {dict_to_table_html(pymses_dict)}
        #     </div>
        #     """

        #     # Return sink_rows + updated states
        #     return status_text, hydro_html, sink_html, pymsesrc_html, sink_rows, cfg, pipe