from radprocess.utils.ramsesinfo import SinkInfo

def get_snapshot_number(input_string):
    """Splits a string at the underscore and returns the part after it."""
    return input_string.split('_')[-1].split('/')[0] if '_' in input_string else None


def hydro_file_descriptor(path):
    """Parses the hydro file descriptor to extract the number of variables, a dictionary of variable names, and the number of dust ratios."""
    filename = path + "hydro_file_descriptor.txt"
    variables = {}
    nvar = 0


    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("nvar"):
                nvar = int(line.split('=')[1].strip())
            elif line.startswith("variable"):
                parts = line.split(":")
                index = int(parts[0].split("#")[1].strip())
                name = parts[1].strip()
                variables[index] = name

    nb_dust_ratios = sum(1 for name in variables.values() if "dust_ratio" in name) #count the number of dust ratios

    return nvar, variables, nb_dust_ratios


def sink_info(path):
    """Parses the sink file to extract the number of sinks and a dictionary of column names and values."""
    
    model_nb = get_snapshot_number(path)
    filename = "sink_" + model_nb + ".info"

    with open(path+filename, 'r') as file:
        lines = file.readlines()

    # Extract the number of sinks from line 1
    num_sinks = int(lines[0].split('=')[1].strip())

    # Extract column names from line 3
    column_names = lines[2].split()

    # Extract values from line 5
    values = lines[4].split()

    # Create a dictionary mapping column names to their respective values
    sink_data = {column_names[i]: values[i] for i in range(len(column_names))}

    #return num_sinks, column_names, sink_data

    return SinkInfo(
        columns=column_names,
        rows=[sink_data],  # list of one row
        num_sinks=num_sinks
    )

