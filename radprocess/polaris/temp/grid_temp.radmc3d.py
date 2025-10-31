# Ilseung Han
#   22.05.2025
#   28.07.2025

import os
import struct
import numpy as np

with open('dust_temperature.bdat', 'rb') as f:
    f.seek(4 * 8)   # skip 4-header entires (each 8 bytes)
    trad = np.fromfile(f, dtype = np.float32)

output_file = 'grid_temp.radmc3d.dat'
os.system('rm %s' % output_file)
with open(output_file, 'wb') as out_f:
    with open('grid_temp.polaris.dat','rb') as f:
        # grid ID
        grid_id = struct.unpack('H', f.read(2))[0]
        out_f.write(struct.pack('H', grid_id))

        # number of parameters in each cell
        n_param = struct.unpack('H', f.read(2))[0]
        out_f.write(struct.pack('H', n_param))

        # parameters ID
        for i in range(n_param):
            param_id = struct.unpack('H', f.read(2))[0]
            out_f.write(struct.pack('H', param_id))
        
        # grid size
        grid_size = struct.unpack('d', f.read(8))[0]
        out_f.write(struct.pack('d', grid_size))
        
        # grid
        leaf_index = 0
        # for i in range(20):
        for i in range(2000000):
            # branch (0) / leaf (1): POLARIS
            # branch (1) / leaf (0): RADMC-3D
            buf = f.read(2)
            if len(buf) < 2:
                print(f"End of file reached at index {i} (Branch/Leaf)")
                break
            t = struct.unpack('H', buf)[0]
            out_f.write(struct.pack('H', t))

            # level
            buf = f.read(2)
            if len(buf) < 2:
                print(f"End of file reached at index {i} (Level)")
                break
            level = struct.unpack('H', buf)[0]
            out_f.write(struct.pack('H', level))

            if t == 1: # only for leaf cells
                # parameters values
                for j in range(n_param):
                    buf = f.read(4)
                    if len(buf) < 4:
                        print(f"End of file reached at grid {i}, param {j}")
                        break
                    if j == 1: # dust temperature
                        val = trad[leaf_index]
                    else:
                        val = struct.unpack('f', buf)[0]
                    out_f.write(struct.pack('f', val))
                leaf_index += 1
print('Results saved in %s' % output_file)
