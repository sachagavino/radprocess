# Ilseung Han
#   16.07.2025
#   18.07.2025
#   22.07.2025
#   26.07.2025

import os

def dustopac(ndust, model):
    os.system('rm -rf dustopac.inp')
    with open('dustopac.inp', 'w') as f:
        f.write('2\n')                          # format number (always 2!)
        f.write('%d' % ndust)                   # the number of dust species
        for i in range(ndust):
            f.write('\n-')
            f.write('\n1')                      # input file style (1: dustkappa_*.inp)
            f.write('\n0')                      # normal thermal grains (always 0!)
            if type(model) == str:
                f.write('\n%s' % model)         # the name of dust species
            else: f.write('\n%s' % model[i])

model = [f'dust_mixture_{(i + 1):03d}' for i in range(10)]
dustopac(ndust = len(model), model = model)
