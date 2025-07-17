# Ilseung Han
#   15.07.2025
#   16.07.2025
#       "This file is only needed if you want to create a spectrum at a special set of wavelengths (otherwise use radmc3d sed)."

import os
import astropy.units as u

def camera_wavelength_micron(wave = [900, 1300, 3100]):
    os.system('rm -rf camera_wavelength_micron.inp')
    with open('camera_wavelength_micron.inp', 'w') as f:
        f.write('%d' % len(wave))
        for w in wave:
            f.write('\n%e' % w)
