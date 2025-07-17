# Ilseung Han
#   26.03.2025
#   16.07.2025

import os
import numpy as np

def wavelength_micron(wave = [0.27, 4000], nwave = 200):
    wave = np.logspace(np.log10(wave[0]), np.log10(wave[1]), nwave)
    os.system('rm -rf wavelength_micron.inp')
    with open('wavelength_micron.inp', 'w') as f:
        f.write('%d' % len(wave))
        for w in wave:
            f.write('\n%e' % w)
