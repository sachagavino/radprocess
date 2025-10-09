# Ilseung Han
#   21.07.2025
#   30.07.2025
#   31.07.2025

import os
import numpy as np

# os.system('python output_00013.star.py')
os.system('rm -rf output_00013.temp.cmd')
with open('output_00013.temp.cmd', 'w') as f:
    # <common>
    f.write('<common>\n')

    # dust components
    size = np.logspace(np.log10(5e-09), np.log10(1e-02), 11)
    path = '/data/pebbles/software/polaris/github/input/dust_nk/silicate_d03.nk'
    # dist = 'plaw'
    # frac = 1.0
    # mden = 0
    for i in range(10):
        f.write('\n\t<dust_component id = "%d"> "%s" "plaw" 1.0 0 %.2e %.2e -3.5'
                % (i, path, size[i], size[i + 1]))
        # WARNING: Number of scattering angles might be too low (max rel diff = 0.901011).
        #   If required, increase 'NANG' or decrease 'MAX_MIE_SCA_REL_DIFF' (for x < 100) in src/Typedefs.h and recompile!

    # thread numbers
    f.write('\n\n\t<nr_threads> -1\n')

    f.write('\n</common>\n')

    # <task> 1
    f.write('\n<task> 1\n')

    # Monte-Carlo calculation
    f.write('\n\t<cmd> CMD_TEMP\n')

    # radiation sources
        # xpos, ypos, zpos = 3.3881627500856327e10, 1.2303378480752747e11, -7.259610992590149e10
        # Note that these values are from "output_00080," not "output_00013."
    xpos, ypos, zpos = -159024056033596.66, 76243296132403.3, -165104488761196.7
    f.write('\n\t<source_star nr_photons = "1"> %e %e %e 1 5804' % (xpos, ypos, zpos))
    f.write('\n\t<source_isrf nr_photons = "0">\n')

    # paths for input/output
    f.write('\n\t<path_grid> "./ramses_grid_00013.dat"')
    f.write('\n\t<path_out> "./"\n')
    # f.write('\n\t<path_out> "./grid_temp.polaris.dat"\n')
    
    # mean molecular weight
    f.write('\n\t<mu> 2.31\n')

    # dust-to-gas mass ratio
    f.write('\n\t<mass_fraction> 0\n')

    f.write('\n</task>')

os.system('rm -rf data/ plots/ grid_temp.dat')
os.system('polaris output_00013.temp.cmd')
os.system('mv grid_temp.dat grid_temp.polaris.dat')
