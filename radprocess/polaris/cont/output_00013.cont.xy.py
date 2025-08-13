# Ilseung Han
#   21.07.2025
#   12.08.2025

import os
import numpy as np
import astropy.units as u

os.system('rm -rf output_00013.cont.xy.cmd')
with open('output_00013.cont.xy.cmd', 'w') as f:
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

    # ray-tracing
    f.write('\n\t<cmd> CMD_DUST_EMISSION\n')

    # paths for input/output
    # f.write('\n\t<path_grid> "./grid_temp.dat"')
    f.write('\n\t<path_grid> "../temp/grid_temp.polaris.dat"')
    f.write('\n\t<path_out> "./xy/"\n')
    
    # mean molecular weight
    f.write('\n\t<mu> 2.31\n')

    # dust-to-gas mass ratio
    f.write('\n\t<mass_fraction> 0\n')

    # alignment mechanism of non-spherical grains
    f.write('\n\t<align> ALIG_PA\n')

    # # rotation axes for the first/second rotation angles
    f.write('\n\t<axis1> 1 0 0\n')  # x-axis
    f.write('<axis2> 0 1 0\n')      # y-axis

    # detector for dust emission
    npix = 301
    wave = 0.9 * u.mm
    dist = 250 * u.pc
    f.write('\n\t<detector_dust nr_pixel = "%d"> %e %e 1 1 0 0 %e\n'
            % (npix, wave.to(u.m).value, wave.to(u.m).value, dist.to(u.m).value))

    # # sub-pixeling
    # f.write('\n\t<max_subpixel_lvl> 3\n')
    # # WARNING: level of subpixeling (3) might be too low
    # #   if required, increase the maximum level of subpixeling with <max_subpixel_lvl> in the command file

    # write midplane temperatures
    f.write('\n\t<write_inp_midplanes> %d' % npix)
    f.write('\n\t<write_3d_midplanes> 1 %d\n' % npix)

    f.write('\n</task>')

os.system('rm -rf xy/')
os.system('polaris output_00013.cont.xy.cmd')
