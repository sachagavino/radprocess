import os

import numpy as np

import radprocess.ramses.read as read


def pymsesrc(ramses_dir):
    """
    Automatically generate ~/.pymses/pymsesrc based on the RAMSES
    file_descriptor.txt found in the simulation output.

    Handles:
      - scalar fields (density, pressure, temperature, dust_ratio_N, fluid_density_N ...)
      - vector fields (velocity_x/y/z, B_left_x/y/z, fluid_v_x/y/z_*, ...)
      - gravity fields (file_type="grav")
    """

    # ----------------------------------------------------------
    # Load descriptors
    # ----------------------------------------------------------
    nvar, variables, nb_ratios = read.hydro_file_descriptor(ramses_dir)
    nvarflui, fluid_variables, nb_fluids = read.other_file_descriptor(ramses_dir)

    # ----------------------------------------------------------
    # Build combined varmap (0-based for PyMSES)
    # ----------------------------------------------------------
    # Hydro first
    varmap = {i - 1: name for i, name in variables.items()}
    offset = len(varmap)

    # Append fluid variables (if any)
    if fluid_variables is not None:
        for name in fluid_variables.values():
            varmap[offset] = name
            offset += 1


    # ----------------------------------------------------------
    # Prepare grouping of vector fields
    # ----------------------------------------------------------
    scalars = {}        # name -> ivar
    vector_groups = {} # basename -> [ivar_x, ivar_y, ivar_z]

    for ivar, name in varmap.items():
        if name.endswith(("_x", "_y", "_z")):
            base = name[:-2]
            axis = name[-1]

            if base not in vector_groups:
                vector_groups[base] = [None, None, None]

            idx = {"x": 0, "y": 1, "z": 2}[axis]
            vector_groups[base][idx] = ivar
        else:
            scalars[name] = ivar

    # ----------------------------------------------------------
    # Remove vector components from scalar list
    # ----------------------------------------------------------
    for base in vector_groups:
        for suffix in ["_x", "_y", "_z"]:
            scalars.pop(base + suffix, None)

    # ------------------ Always include density first -----------------
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

        # ----------------------- Vector fields -----------------------
        for base, ivars in vector_groups.items():
            if None not in ivars:
                f.write(",\n")
                f.write(
                    f'            {{"__type__": "vector_field", "__file_type__": "hydro", '
                    f'"name": "{base}", "ivars": {ivars}}}'
                )

        # ----------------------- Scalar fields -----------------------
        for name, ivar in scalars.items():
            f.write(",\n")

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
