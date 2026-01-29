import os
import re
from collections import defaultdict

import radprocess.ramses.read as read


def pymsesrc(ramses_dir):

    # ----------------------------------------------------------
    # Regex
    # ----------------------------------------------------------
    # hydro vectors: velocity_x, B_left_y, etc.
    std_vec_re = re.compile(r"^(.*)_(x|y|z)$")

    # fluid vectors: fluid_v_x_17
    fluid_vec_re = re.compile(r"^(fluid_v)_(x|y|z)_(\d+)$")

    # ----------------------------------------------------------
    # Load descriptors
    # ----------------------------------------------------------
    _, hydro_vars, _ = read.hydro_file_descriptor(ramses_dir)
    _, fluid_vars, _ = read.other_file_descriptor(ramses_dir)

    # Convert to 0-based local indexing
    hydro_map = {i - 1: name for i, name in hydro_vars.items()}
    fluid_map = {i - 1: name for i, name in fluid_vars.items()} if fluid_vars else {}

    print("Hydro ivars:", min(hydro_map), max(hydro_map))
    if fluid_map:
        print("Fluid ivars:", min(fluid_map), max(fluid_map))

    # ----------------------------------------------------------
    # Containers
    # ----------------------------------------------------------
    hydro_scalars = {}
    hydro_vectors = defaultdict(lambda: [None, None, None])

    fluid_scalars = {}
    fluid_vectors = defaultdict(lambda: [None, None, None])

    # ----------------------------------------------------------
    # Parse hydro variables
    # ----------------------------------------------------------
    for ivar, name in hydro_map.items():

        m = std_vec_re.match(name)
        if m:
            base, axis = m.groups()
            hydro_vectors[base][{"x":0,"y":1,"z":2}[axis]] = ivar
        else:
            hydro_scalars[name] = ivar

    # ----------------------------------------------------------
    # Parse fluid variables
    # ----------------------------------------------------------
    for ivar, name in fluid_map.items():

        m = fluid_vec_re.match(name)
        if m:
            base, axis, idx = m.groups()
            key = f"{base}_{idx}"           # fluid_v_17
            fluid_vectors[key][{"x":0,"y":1,"z":2}[axis]] = ivar
        else:
            fluid_scalars[name] = ivar

    # ----------------------------------------------------------
    # Remove vector components from scalar dicts
    # ----------------------------------------------------------
    for base, ivars in hydro_vectors.items():
        for iv in ivars:
            if iv is not None:
                hydro_scalars.pop(hydro_map[iv], None)

    for base, ivars in fluid_vectors.items():
        for iv in ivars:
            if iv is not None:
                fluid_scalars.pop(fluid_map[iv], None)

    # ----------------------------------------------------------
    # Density first (hydro)
    # ----------------------------------------------------------
    density_ivar = hydro_scalars.pop("density", 0)

    # ----------------------------------------------------------
    # Write ~/.pymses/pymsesrc
    # ----------------------------------------------------------
    pymses_dir = os.path.expanduser("~/.pymses")
    os.makedirs(pymses_dir, exist_ok=True)

    rc_file = os.path.join(pymses_dir, "pymsesrc")
    print(f"Writing pymsesrc → {rc_file}")

    with open(rc_file, "w") as f:

        f.write('{\n')
        f.write(' "Version": 1,\n')
        f.write(' "Multiprocessing max. nproc": 8,\n')
        f.write(' "RAMSES": {\n')
        f.write('  "ndimensions": 3,\n')
        f.write('  "amr_field_descr": [\n')

        # ---------------- Density ----------------
        f.write(
            f'   {{"__type__":"scalar_field","__file_type__":"hydro","name":"density","ivar":{density_ivar}}}'
        )

        # ---------------- Hydro vectors ----------------
        for base, ivars in hydro_vectors.items():
            if None not in ivars:
                f.write(",\n")
                f.write(
                    f'   {{"__type__":"vector_field","__file_type__":"hydro","name":"{base}","ivars":{ivars}}}'
                )

        # ---------------- Hydro scalars ----------------
        for name, ivar in hydro_scalars.items():
            f.write(",\n")
            f.write(
                f'   {{"__type__":"scalar_field","__file_type__":"hydro","name":"{name}","ivar":{ivar}}}'
            )

        # ---------------- Fluid vectors ----------------
        for base, ivars in fluid_vectors.items():
            if None not in ivars:
                f.write(",\n")
                f.write(
                    f'   {{"__type__":"vector_field","__file_type__":"mf","name":"{base}","ivars":{ivars}}}'
                )

        # ---------------- Fluid scalars ----------------
        for name, ivar in fluid_scalars.items():
            f.write(",\n")
            f.write(
                f'   {{"__type__":"scalar_field","__file_type__":"mf","name":"{name}","ivar":{ivar}}}'
            )

        f.write("\n  ]\n }\n}\n")

    return rc_file

