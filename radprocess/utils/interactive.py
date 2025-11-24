import gradio as gr
from dataclasses import fields, is_dataclass

def interactive_config_gradio(dataclass_obj):
    """
    Create a Gradio interface for a dataclass with sections INOUT, PYMSESRC, SIM.
    Booleans are checkboxes, strings/numbers are editable fields. Updates dataclass in-place.
    """
    if not is_dataclass(dataclass_obj):
        raise TypeError("Expected a dataclass instance.")

    # Prepare a dictionary mapping tab names to lists of input components
    tabs_content = {}

    for f in fields(dataclass_obj):
        val = getattr(dataclass_obj, f.name)
        title = f.name.upper()

        if is_dataclass(val):
            inputs = []
            outputs = []

            for subf in fields(val):
                subval = getattr(val, subf.name)
                desc = subf.metadata.get("desc", "")
                typ = subf.type

                # Boolean → Checkbox
                if typ is bool:
                    inp = gr.Checkbox(value=subval, label=subf.name)
                # Integer → Number
                elif typ is int:
                    inp = gr.Number(value=subval, label=subf.name, precision=0)
                # Float → Number
                elif typ is float:
                    inp = gr.Number(value=subval, label=subf.name)
                # String → Textbox
                elif typ is str:
                    width = 700 if "dir" in subf.name.lower() else 200
                    inp = gr.Textbox(value=subval, label=subf.name, lines=1)
                # Fallback
                else:
                    inp = gr.Textbox(value=str(subval), label=subf.name)

                inputs.append(inp)
                outputs.append(subf.name)

            tabs_content[title] = inputs

    # Function to update the dataclass from Gradio inputs
    def update_dataclass(*args):
        idx = 0
        for f in fields(dataclass_obj):
            val = getattr(dataclass_obj, f.name)
            if is_dataclass(val):
                for subf in fields(val):
                    setattr(val, subf.name, args[idx])
                    idx += 1
        return "Configuration updated!"

    # Flatten all inputs to feed Gradio function
    all_inputs = sum(tabs_content.values(), [])

    with gr.Blocks() as demo:
        gr.Markdown("### Editable Configuration")
        tab_names = list(tabs_content.keys())
        for tname in tab_names:
            with gr.Tab(tname):
                for comp in tabs_content[tname]:
                    comp.render()
        # Submit button
        submit_btn = gr.Button("Update Configuration")
        output_text = gr.Textbox(value="", label="Status", interactive=False)
        submit_btn.click(update_dataclass, all_inputs, output_text)

    demo.launch(share=True)


if __name__ == "__main__":
    # Example usage
    from config import ConfigParams  # your dataclass
    cfg = ConfigParams()
    interactive_config_gradio(cfg)
