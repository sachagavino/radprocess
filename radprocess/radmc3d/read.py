import numpy as np


def grid(path, filename=None):

    if (filename == None):
        filename = path + "amr_grid.inp"
 


    f = open(filename,"r")


    f.readline()
    grid_type = int(f.readline())

    if grid_type == 0:
        coordsystem = int(f.readline())
        f.readline()
        incl_x, incl_y, incl_z = [int(x) for x in f.readline().split()]
        nx, ny, nz = [int(x) for x in f.readline().split()]

        grid = []

        # x, y, z = [], [], []

        # for ix in range(nx+1):
        #     x.append(float(f.readline()))
        # for iy in range(ny+1):
        #     y.append(float(f.readline()))
        # for iz in range(nz+1):
        #     z.append(float(f.readline()))

        x = np.array(f.readline().replace("\n","").split(), dtype=float)
        y = np.array(f.readline().replace("\n","").split(), dtype=float)
        z = np.array(f.readline().replace("\n","").split(), dtype=float)

    else:
        print('grid type must be regular for the moment.')
        exit(1)

    f.close()

    return nx, ny, nz, x, y, z


def stars(path, filename=None):
    """
    WARNING: not working more than one star and not working for binary files yet.
    """

    if (filename == None):
        filename = path + "stars.inp"
 
    f = open(filename,"r")

    f.readline()
    nstars, nb_lam = [int(x) for x in f.readline().split()]
    r_star, m_star, x_star, y_star, z_star = [float(x) for x in f.readline().split()]

    wave = np.empty((nb_lam,))
    spectrum = np.empty((nb_lam,))

    for j in range(0, nb_lam, 1):
        wave[j] = float(f.readline())

    first_line = float(f.readline())

    if first_line < 0.:
        Tstar = -first_line
        spectrum=-1
    else:
        Tstar=-1
        spectrum[0]=first_line
        for j in range(1, nb_lam, 1):
            spectrum[j] = float(f.readline())

    return nb_lam, wave, r_star, m_star, Tstar, spectrum




def dust_density(path, filename=None, ext=None, binary=False):

    if (filename == None):
        if (ext == None):
            if binary:
                filename = path + "dust_density.binp"
            else:
                filename = path + "dust_density.inp"
        else:
            if binary:
                filename = path + "dust_density_"+str(ext)+".binp"
            else:
                filename = path + "dust_density_"+str(ext)+".inp"

    if binary:
        f = open(filename, "rb")
        data = np.fromfile(filename)
    else:
        f = open(filename,"r")

    if binary:
        int.from_bytes(f.read(8), byteorder="little")
        int.from_bytes(f.read(8), byteorder="little")
        ncells = int.from_bytes(f.read(8), byteorder="little")
        nspecies = int.from_bytes(f.read(8), byteorder="little")
    else:
        f.readline()
        ncells = int(f.readline())
        nspecies = int(f.readline())

    density = []
    index = 3
    #for i in range(nspecies):
    dens = np.empty((nspecies*ncells,))

    for j in range(nspecies*ncells):
        if binary:
            dens[j] = data[index+j]
        else:
            dens[j] = float(f.readline())

    density.append(dens)

    f.close()

    return density

def dust_temperature(path, filename=None, ext=None, binary=False):

    if (filename == None):
        if (ext == None):
            if binary:
                filename = path + "dust_temperature.bdat"
            else:
                filename = path + "dust_temperature.dat"
        else:
            if binary:
                filename = path + "dust_temperature_"+str(ext)+".bdat"
            else:
                filename = path + "dust_temperature_"+str(ext)+".dat"

    if binary:
        try:
            f = open(filename, "rb")
            data = np.fromfile(filename)
        except IOError:
            return []

    else:
        try:
            f = open(filename,"r")
        except IOError:
            return []


    if binary:
        int.from_bytes(f.read(8), byteorder="little")
        int.from_bytes(f.read(8), byteorder="little")
        ncells = int.from_bytes(f.read(8), byteorder="little")
        nspecies = int.from_bytes(f.read(8), byteorder="little")
    else:
        f.readline()
        ncells = int(f.readline())
        nspecies = int(f.readline())

    temperature = []
    index = 3
    #for i in range(nspecies):
    temp = np.empty((nspecies*ncells,))

    for j in range(nspecies*ncells):
        if binary:
            temp[j] = data[index+j]
        else:
            temp[j] = float(f.readline())

    temperature.append(temp)

    f.close()

    return temperature


def localfield(path, filename=None):

    # if (filename == None):
    #     filename = path + "mean_intensity.out"


    # f = open(filename,"r")

    # f.readline()
    # ncells = int(f.readline())
    # nfreq = int(f.readline())
    # f.readline()

    # localfield = []
    # # for i in range(nfreq*ncells):
    # #     localfield.append(float(f.readline()))

    # for i in range(nfreq):
    #     field = np.empty((ncells,))

    #     for j in range(ncells):
    #         field[j] = float(f.readline())

    #     localfield.append(field)

    # f.close()

    # return localfield


    if (filename == None):
        filename = path + "mean_intensity.out"

    f = open(filename,"r")
  
    f.readline()
    nb_pt = int(f.readline())
    nlam_mono = int(f.readline())
    wave_mono = np.empty((nlam_mono,))
    localfield = np.empty((nlam_mono*nb_pt,))

    #for nbl in range(0, nlam_mono):
    wave_mono = [float(x) for x in f.readline().split()] #in Hz
    wave_mono = np.asarray(wave_mono)

    for j in range(0, nlam_mono*nb_pt, 1):
        localfield[j] = float(f.readline())
  
    return nlam_mono, wave_mono, localfield


def mcmono(path, filename=None):

    if (filename == None):
        filename = path + "mcmono_wavelength_micron.inp"
    
    f = open(filename,"r")

    nblam_mono = int(f.readline())
    wave_mono = np.empty((nblam_mono,))

    for nbl in range(0, nblam_mono):
        wave_mono[nbl] = f.readline().split()
  
    return nblam_mono, wave_mono



def external_source(path, filename=None):
    """
    WARNING: not working more than one star and not working for binary files yet.
    """

    if (filename == None):
        filename = path + "external_source.inp"
 
    try:
        f = open(filename,"r")

        f.readline()
        nb_lam = int(f.readline())

        wave = np.empty((nb_lam,))
        spectrum = np.empty((nb_lam,))

        for j in range(0, nb_lam, 1):
            wave[j] = float(f.readline())

        for j in range(0, nb_lam, 1):
            spectrum[j] = float(f.readline())
    except FileNotFoundError:
        spectrum = []

    return spectrum
