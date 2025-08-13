# Ilseung Han
#   13.05.2025
#   28.07.2025

import os
import struct
import numpy as np

output_file = 'dust_temperature.bdat.txt'
os.system('rm -rf %s' % output_file)
with open(output_file, 'w') as out_f:
    with open('dust_temperature.bdat','rb') as f:
        # iformat
        iformat = struct.unpack('q', f.read(8))[0]
        print(iformat, file = out_f)

        # precision
        precis = struct.unpack('q', f.read(8))[0]
        print(precis, file = out_f)

        # number of cells
        nrcells = struct.unpack('q', f.read(8))[0]
        print(nrcells, file = out_f)

        # number of dust species
        nrspec = struct.unpack('q', f.read(8))[0]
        print(nrspec, file = out_f)

        # temperature
        # for i in range(100):
        for i in range(nrcells * nrspec):
            temp = struct.unpack('f', f.read(4))[0]
            # 4-byte temperature in POLARIS
            print(temp, file = out_f)
print('Results saved in %s' % output_file)

# temp = np.loadtxt(output_file, skiprows = 4)
# temp = temp.reshape([-1, nrspec], order = 'F')
