# Ilseung Han
#   18.07.2025
#   22.07.2025
#   30.07.2025

import os
import numpy as np

for i in range(10):
    os.system('rm -rf data/dustkappa_dust_mixture_%03d.inp' % (i + 1))
    with open('dustkappa_dust_mixture_%03d.inp' % (i + 1), 'w') as f:
        dust = np.loadtxt('data/dust_mixture_%03d.dat' % (i + 1), skiprows = 27)
        wave = dust[:,  0] / 1e-06  # wavelength in micrometers
        kabs = dust[:, -4] * 1e+01  # absorption opacity
        ksca = dust[:, -2] * 1e+01  # scattering opacity
        f.write('2\n')              # format number (2: isotropic scattering)
        f.write('%d' % len(wave))   # the number of wavelengths
        for j in range(len(wave)):
            f.write('\n%e %e %e' % (wave[j], kabs[j], ksca[j]))
