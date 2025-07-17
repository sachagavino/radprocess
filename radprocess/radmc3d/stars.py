# Ilseung Han
#   26.03.2025
#   16.07.2025

import os
import numpy as np

def stars(nstar, pstar, rstar, mstar, tstar, iformat = 2, wave = [0.27, 4000], nwave = 200):
    wave = np.logspace(np.log10(wave[0]), np.log10(wave[1]), nwave)
    os.system('rm -rf stars.inp')
    with open('stars.inp', 'w') as f:
        f.write('%d\n' % iformat)           # format number
                                            #   1: frequencies; 2: wavelengths
        f.write('%d '  % nstar)             # the number of stars
        f.write('%d\n' % len(wave))         # the number of wavelengths
        for i in range(nstar):
            f.write('%e '  % rstar[i])      # stellar radius (cm)
            f.write('%e '  % mstar[i])      # stellar mass (g)
            f.write('%e '  % pstar[i][0])   # stellar position X (cm)
            f.write('%e '  % pstar[i][1])   # stellar position Y (cm)
            f.write('%e\n' % pstar[i][2])   # stellar position Z (cm)
        for i in range(nstar):
            for w in wave:
                f.write('\n%e' % w)         # emitting wavelength
        for t in tstar:
            f.write('\n-%d' % t)            # stellar temperature
