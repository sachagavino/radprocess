# Ilseung Han
#   16.07.2025

import os

def radmc3d(nphot = 5000000):
    os.system('rm -rf radmc3d.inp')
    with open('radmc3d.inp', 'w') as f:
        f.write('nphot = %d\n' % nphot)         # the number of photons for thermal MC (for temperature)
        f.write('nphot_scat = 1000000\n')       # the number of photons for scattering MC (for imaging)
        f.write('scattering_mode = 1\n')        # isotropic scattering
        f.write('scattering_mode_max = 1\n')    # isotropic scattering
        f.write('setthreads = 4\n')             # the number of threads for MC calculations
        f.write('modified_random_walk = 1')     # modified random walk (MRW)
