import os

import numpy as np


def pymsesrc(self, 
                    ndust=1,
                    rho=True, #always true
                    dustratios=True, 
                    vel=True, 
                    bl=False, 
                    br=False, 
                    p=False, 
                    xi=False, 
                    phi=False, 
                    g=False):
    # Initialize the variable i
    i = 0
    # Define the pymsesrc path
    pymses_directory = os.path.expanduser("~/.pymses")
    
    # Check if the directory exists, and create it if it doesn't
    if not os.path.exists(pymses_directory):
        os.makedirs(pymses_directory)

    # Define the file path in the pymses directory
    rc_file = os.path.join(pymses_directory, "pymsesrc")
    print(rc_file)

    # Create and write to the file
    with open(rc_file, "w") as f:
        f.write('{\n')
        f.write('    "Version": 1,\n')
        f.write('    "Multiprocessing max. nproc": 8,\n')
        f.write('    "RAMSES":{\n')
        f.write('        "ndimensions": 3,\n')
        f.write('        "amr_field_descr": [\n')
        f.write('            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "rho", "ivar": 0}')
        if (dustratios==True):
            while i < ndust:
                f.write(',\n')
                f.write('            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "dustratio%d", "ivar": %d}' % (i+1, i+11))
                i += 1
        if (vel==True):
            f.write(',\n')
            f.write('            {"__type__": "vector_field", "__file_type__": "hydro", "name": "vel", "ivars": [1, 2, 3]}')
        if (bl==True):
            f.write(',\n')
            f.write('            {"__type__": "vector_field", "__file_type__": "hydro", "name": "Bl", "ivars": [4, 5, 6]}')
        if (br==True):
            f.write(',\n')
            f.write('            {"__type__": "vector_field", "__file_type__": "hydro", "name": "Br", "ivars": [7, 8, 9]}')
        if (p==True):
            f.write(',\n')
            f.write('            {"__type__": "scalar_field", "__file_type__": "hydro", "name": "P", "ivar": 10}')
        if (xi==True):
            f.write(',\n')
            f.write('            {"__type__": "multivalued_field", "__file_type__": "hydro", "name": "Xi", "ivar_first": 11, "nb_vars": 7}')
        if (phi==True):
            f.write(',\n')
            f.write('            {"__type__": "scalar_field", "__file_type__": "grav", "name": "phi", "ivar": 0}')
        if (g==True):
            f.write(',\n')
            f.write('            {"__type__": "vector_field", "__file_type__": "grav", "name": "g", "ivars": [1, 2, 3]}')
        f.write('\n')
        f.write('        ]\n')
        f.write('    }\n')
        f.write('}\n')