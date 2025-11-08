import os

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

