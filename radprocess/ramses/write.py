import os

import numpy as np

import radprocess.ramses.read as read


def pymsesrc(ramses_dir):
    """
    Automatically generate ~/.pymses/pymsesrc based on the RAMSES
    hydro_file_descriptor.txt found in the simulation output.

    Handles:
      - scalar fields (density, pressure, temperature, dust_ratio_N, ...)
      - vector fields (velocity_x/y/z, B_left_x/y/z, ...)
      - gravity fields (file_type="grav")
    """

    # ----------------------------------------------------------
    # Load descriptor
    # ----------------------------------------------------------
    nvar, variables, nb_dust = read.hydro_file_descriptor(ramses_dir)

    # variables: {1: "density", 2: "velocity_x", ...}
    # Convert to 0-based indexing for PyMSES
    varmap = {i - 1: name for i, name in variables.items()}

    # ----------------------------------------------------------
    # Prepare grouping of vector fields
    # ----------------------------------------------------------
    scalars = {}           # name -> ivar
    vector_groups = {}     # basename -> [ivar_x, ivar_y, ivar_z]

    for ivar, name in varmap.items():
        # Identify vector components: something ending with _x, _y, or _z
        if name.endswith(("_x", "_y", "_z")):
            base = name[:-2]  # remove _x, _y, _z
            axis = name[-1]   # x, y, z

            if base not in vector_groups:
                vector_groups[base] = [None, None, None]

            idx = {"x": 0, "y": 1, "z": 2}[axis]
            vector_groups[base][idx] = ivar

        else:
            scalars[name] = ivar

    # ----------------------------------------------------------
    # Remove parts of vector fields from scalars dict
    # ----------------------------------------------------------
    for base in vector_groups:
        for suffix in ["_x", "_y", "_z"]:
            full = base + suffix
            if full in scalars:
                del scalars[full]

    # ------------------ Always include density -----------------
    # density is always var #1 → ivar 0
    # But if the descriptor included it already, ensure it's first.
    if "density" in scalars:
        density_ivar = scalars.pop("density")
    else:
        density_ivar = 0

    # ----------------------------------------------------------
    # Path to ~/.pymses/pymsesrc
    # ----------------------------------------------------------
    pymses_dir = os.path.expanduser("~/.pymses")
    os.makedirs(pymses_dir, exist_ok=True)

    rc_file = os.path.join(pymses_dir, "pymsesrc")
    print(f"Writing pymsesrc → {rc_file}")

    # ----------------------------------------------------------
    # Write JSON manually (PyMSES legacy format)
    # ----------------------------------------------------------
    with open(rc_file, "w") as f:
        f.write('{\n')
        f.write('    "Version": 1,\n')
        f.write('    "Multiprocessing max. nproc": 8,\n')
        f.write('    "RAMSES": {\n')
        f.write('        "ndimensions": 3,\n')
        f.write('        "amr_field_descr": [\n')

        # ----------------------- Density first -----------------------
        f.write(
            f'            {{"__type__": "scalar_field", '
            f'"__file_type__": "hydro", "name": "density", "ivar": {density_ivar}}}'
        )

        first_written = True

        # ----------------------- Vector fields -----------------------
        for base, ivars in vector_groups.items():
            if None not in ivars:  # only write if full vector exists
                f.write(",\n")
                f.write(
                    f'            {{"__type__": "vector_field", "__file_type__": "hydro", '
                    f'"name": "{base}", "ivars": {ivars}}}'
                )

        # ----------------------- Scalar fields -----------------------
        for name, ivar in scalars.items():
            f.write(",\n")

            # file type: hydro except gravity
            ftype = "grav" if name.startswith("gravity") else "hydro"

            f.write(
                f'            {{"__type__": "scalar_field", "__file_type__": "{ftype}", '
                f'"name": "{name}", "ivar": {ivar}}}'
            )

        f.write("\n")
        f.write("        ]\n")
        f.write("    }\n")
        f.write("}\n")

    return rc_file