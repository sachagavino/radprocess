# Ilseung Han
#   16.07.2025
#   28.07.2025
#   31.07.2025

import os

def radmc3d(nphot = 1000000):
    os.system('rm -rf radmc3d.inp')
    with open('radmc3d.inp', 'w') as f:
        f.write('nphot = %d\n' % nphot)         # the number of photons for thermal MC (for temperature)
        f.write('nphot_scat = 1000000\n')       # the number of photons for scattering MC (for imaging)
        f.write('scattering_mode = 1\n')        # isotropic scattering
        f.write('scattering_mode_max = 1\n')    # isotropic scattering
        f.write('setthreads = 4\n')             # the number of threads for MC calculations
                                                # Will this parameter have to be a free parameter later?
        f.write('modified_random_walk = 1\n')   # modified random walk (MRW)
        f.write('rto_style = 3\n')              # output in binary form
        f.write('rto_single = 1')               # single precision (same as POLARIS)

radmc3d(nphot = 1000000)
# os.system('radmc3d mctherm')
