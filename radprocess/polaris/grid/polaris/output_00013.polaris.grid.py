# Ilseung Han
#   18.07.2025
#   30.07.2025

from ram2pol import *

datapath    = '/nfs/pic.es/user/i/ihan/disk/anaelle/ENYGMA/RAMSES_MODELS/TESTMODEL_NODUSTEVOL'
num         = '13'
outpath     = './'
starpos     = [-159024056033596.66, 76243296132403.3, -165104488761196.7]
    # output_00013.starpos.py
# convert_ramses2polaris(datapath, num, outpath, starpos)
convert_ramses2polaris(datapath, num, outpath, starpos, size_hole_au = 8)
