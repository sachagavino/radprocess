import glob
import os
import re
#import json

import numpy as np

from radprocess.utils.ramsesinfo import SinkInfo

def get_snapshot_number(input_string):
    """Splits a string at the underscore and returns the part after it."""
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



def sink_info(path):
    """Parses the sink file to extract information about sinks."""
    
    model_nb = get_snapshot_number(path)
    filename = f"sink_{model_nb}.info"
    filename2 = f"sink_{model_nb}.csv"
    filepath = path + filename
    filepath2 = path + filename2

    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Extract the number of sinks from line 1
    num_sinks = int(lines[0].split('=')[1].strip())

    # --- Column names on line 3 (index 2) ---
    column_sinks = lines[2].split()

    sinks = np.loadtxt(filepath2, delimiter=',')

    rows = []

    # --- Parse all sink entries ---
    # sink lines start at line 4 (index 4) and end before last separator line
    for line in lines[4:]:
        line = line.strip()

        # stop before the last line of ====== separator
        if line.startswith("===="):
            break

        if not line:
            continue

        parts = line.split()
        if len(parts) != len(column_sinks):
            # Some sink lines may include ******** fields → still readable
            # Ensure alignment by trimming or padding
            parts = parts[:len(column_sinks)]

        rows.append({column_sinks[i]: parts[i] for i in range(len(column_sinks))})


    return SinkInfo(
        columns=column_sinks,
        rows=rows,  # list of one row
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
