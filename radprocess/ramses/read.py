import glob
import os
import re
#import json

import numpy as np

from radprocess.utils.ramsesinfo import SinkInfo

def get_snapshot_number(input_string):
    input_string = str(input_string)   # <-- normalize here
    return input_string.split('_')[-1].split('/')[0] if '_' in input_string else None

def hydro_file_descriptor(path):
    """
    Parse the RAMSES hydro_file_descriptor.txt file.

    Returns:
        nvar (int) : total number of hydro variables
        variables (dict): {index: variable_name}
        nb_dust_ratios (int): number of dust ratio variables

    Raises:
        FileNotFoundError: if the file does not exist
    """

    # Ensure path ends with slash
    if not path.endswith("/"):
        path = path + "/"

    filename = path + "hydro_file_descriptor.txt"

    #  Check that the file exists
    if not os.path.isfile(filename):
        raise FileNotFoundError(
            f"hydro_file_descriptor.txt not found in: {path}\n"
            f"Please, check that the RAMSES directory is correct."
        )

    variables = {}
    nvar = 0

    #  Read and parse the file
    try:
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()

                if line.startswith("nvar"):
                    nvar = int(line.split("=")[1].strip())

                elif line.startswith("variable"):
                    parts = line.split(":")
                    index = int(parts[0].split("#")[1].strip())
                    name = parts[1].strip()
                    variables[index] = name

    except Exception as e:
        raise RuntimeError(
            f"Error while reading {filename}:\n{e}\n"
            "Please, check the file hydro_file_descriptor.txt."
        )

    #  Count dust ratios
    nb_dust_ratios = sum("dust_ratio" in name for name in variables.values())

    return nvar, variables, nb_dust_ratios

def other_file_descriptor(path):
    """
    Parse any RAMSES *_file_descriptor.txt except hydro_file_descriptor.txt.

    Returns:
        nvar (int)
        variables (dict): {1-based_index: variable_name}
        nb_fluids (int): inferred number of fluids
    """

    if not path.endswith("/"):
        path += "/"

    # Find any *_file_descriptor.txt except hydro
    candidates = glob.glob(os.path.join(path, "*_file_descriptor.txt"))
    candidates = [f for f in candidates if not f.endswith("hydro_file_descriptor.txt")]

    if not candidates:
        return None, None, 0

    filename = candidates[0]  # assume only one
    print(f"Found extra descriptor: {os.path.basename(filename)}")

    variables = {}
    nvar = 0
    counter = 1  # we assign indices ourselves

    try:
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()

                # nvarflu or similar
                if "=" in line and line.lower().startswith("nvar"):
                    nvar = int(line.split("=")[1].strip())

                elif line.startswith("variable"):
                    # everything after colon is the name
                    name = line.split(":", 1)[1].strip()

                    variables[counter] = name
                    counter += 1

    except Exception as e:
        raise RuntimeError(f"Error reading {filename}:\n{e}")

    # -------------------------------------------------
    # Infer number of fluids safely:
    # count fluid_density_*
    # -------------------------------------------------
    density_vars = [v for v in variables.values() if v.startswith("fluid_density")]

    nb_fluids = len(density_vars)

    return nvar, variables, nb_fluids


def mf_grain_sizes(path):
    """
    Parse the RAMSES multi-fluid info_mf_XXXXX.txt file.

    Returns:
        dict with
            sizes_cm (ndarray): grain size of each fluid (mf_grain_size_i) [cm]
            densities_gcc (ndarray): grain bulk density of each fluid (mf_grain_density_i) [g/cm^3]
            filename (str): path of the parsed file
        or None if no info_mf_*.txt file exists (not a multi-fluid output).
    """

    candidates = sorted(glob.glob(os.path.join(str(path), "info_mf_*.txt")))
    if not candidates:
        return None

    filename = candidates[0]
    entry_re = re.compile(r"^mf_grain_(size|density)_(\d+)\s*=\s*(\S+)")
    values = {"size": {}, "density": {}}

    with open(filename, "r") as f:
        for line in f:
            m = entry_re.match(line.strip())
            if m:
                kind, idx, val = m.groups()
                # Fortran double-precision exponents (1.0d-3)
                values[kind][int(idx)] = float(val.replace("d", "e").replace("D", "E"))

    if not values["size"]:
        raise RuntimeError(f"No mf_grain_size_* entries found in {filename}.")

    order = sorted(values["size"])
    return {
        "sizes_cm": np.array([values["size"][i] for i in order]),
        "densities_gcc": np.array([values["density"].get(i, np.nan) for i in order]),
        "filename": filename,
    }


def namelist_value(path, key):
    """Return the float value of `key` in the RAMSES namelist.txt of `path`, or None."""

    filename = os.path.join(str(path), "namelist.txt")
    if not os.path.isfile(filename):
        return None

    key_re = re.compile(rf"^{re.escape(key)}\s*=\s*([0-9.eEdD+-]+)", re.IGNORECASE)
    with open(filename, "r") as f:
        for line in f:
            m = key_re.match(line.strip())
            if m:
                return float(m.group(1).replace("d", "e").replace("D", "E"))
    return None



from pathlib import Path
import numpy as np

def sink_info(path):
    """Parses the sink file to extract information about sinks."""

    path = Path(path)

    model_nb = get_snapshot_number(path)
    filename  = f"sink_{model_nb}.info"
    filename2 = f"sink_{model_nb}.csv"

    filepath  = path / filename
    filepath2 = path / filename2

    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Extract the number of sinks
    num_sinks = int(lines[0].split('=')[1].strip())

    # Column names
    column_sinks = lines[2].split()

    # Load numeric sink data
    sinks = np.loadtxt(filepath2, delimiter=',')

    # --- FORCE 2D ---
    if sinks.ndim == 1:
        sinks = sinks[np.newaxis, :]

    rows = []

    # Parse sink entries from text file
    for line in lines[4:]:
        line = line.strip()

        if line.startswith("===="):
            break
        if not line:
            continue

        parts = line.split()
        parts = parts[:len(column_sinks)]

        rows.append({column_sinks[i]: parts[i]
                     for i in range(len(column_sinks))})

    # Optional safety check
    assert sinks.shape[0] == num_sinks

    return SinkInfo(
        columns=column_sinks,
        rows=rows,
        num_sinks=num_sinks,
        data=sinks
    )




# def pymsesrc():
#     """Parses the pymsesrc file and returns it as a Python dictionary."""
    
#     pymses_directory = os.path.expanduser("~/.pymses")

#     # Check if the directory exists, and create it if it doesn't
#     if not os.path.exists(pymses_directory):
#         os.makedirs(pymses_directory)

#     filename = os.path.join(pymses_directory, "pymsesrc")

#     # Check that the file exists
#     if not os.path.isfile(filename):
#         raise FileNotFoundError(
#             f"~/.pymses/pymsesrc file not found.\n"
#         )

#     # Read and parse the file
#     try:
#         with open(filename, "r") as f:
#             file_content = json.load(f)  # parse JSON content

#     except json.JSONDecodeError as e:
#         raise RuntimeError(
#             f"Error decoding JSON in {filename}:\n{e}\n"
#             "Please check that the file contains valid JSON."
#         )

#     except Exception as e:
#         raise RuntimeError(
#             f"Error while reading {filename}:\n{e}\n"
#             "Please, check the file pymsesrc."
#         )

#     return file_content
