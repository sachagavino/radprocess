# Ilseung Han
#   26.03.2025
#   16.07.2025
#   22.07.2025
#   31.07.2025

import os
import numpy as np
import astropy.units as u

# def stars(nstar, pstar, rstar, mstar, tstar, iformat = 2, wave = [0.27, 4000], nwave = 200):
def stars(nstar, pstar, rstar, mstar, tstar, iformat = 2, wave = [0.27, 3000], nwave = 200):
    wave = np.logspace(np.log10(wave[0]), np.log10(wave[1]), nwave)
    os.system('rm -rf stars.inp')
    with open('stars.inp', 'w') as f:
        f.write('%d\n' % iformat)       # format number
                                        #   1: frequencies; 2: wavelengths
        f.write('%d '  % nstar)         # the number of stars
        f.write('%d\n' % len(wave))     # the number of wavelengths
        if nstar == 1:
            f.write('%e '  % rstar)     # stellar radius (cm)
            f.write('%e '  % mstar)     # stellar mass (g)
            f.write('%e '  % pstar[0])  # stellar position X (cm)
            f.write('%e '  % pstar[1])  # stellar position Y (cm)
            # f.write('%e\n' % pstar[2])
            f.write('%e' % pstar[2])    # stellar position Z (cm)
            for w in wave:
                f.write('\n%e' % w)     # emitting wavelength
            f.write('\n-%d' % tstar)    # stellar temperature
        else:
            for i in range(nstar):
                f.write('%e '  % rstar[i])
                f.write('%e '  % mstar[i])
                f.write('%e '  % pstar[i][0])
                f.write('%e '  % pstar[i][1])
                # f.write('%e\n' % pstar[i][2])
                f.write('%e' % pstar[i][2])
            for i in range(nstar):
                for w in wave:
                    f.write('\n%e' % w)
            for t in tstar:
                f.write('\n-%d' % t)

nstar = 1

xstar = -159024056033596.66 * u.m
xstar = xstar.to(u.cm).value

ystar = 76243296132403.3 * u.m
ystar = ystar.to(u.cm).value

zstar = -165104488761196.7 * u.m
zstar = zstar.to(u.cm).value

pstar = [xstar, ystar, zstar]

rstar = 1 * u.R_sun
rstar = rstar.to(u.cm).value
    # The radius will have to be the same as the value adopted in POLARIS.

mstar = 1 * u.M_sun
mstar = mstar.to(u.g).value
    # The mass will also have to be the same as the one adopted in POLARIS.

tstar = 5804

stars(nstar, pstar, rstar, mstar, tstar)
