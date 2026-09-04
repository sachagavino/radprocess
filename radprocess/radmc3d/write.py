from pathlib import Path

def stars(radpath,  nstars, masses, pos, radii, lam, teff):


    radpath = Path(radpath)
    filepath = radpath / "stars.inp"
    #nstars = len(rstar)
    nlam = len(lam)

    with open(filepath, "w") as f:

        f.write(str(2)+"\n")
        f.write("{0:d}  {1:d}\n".format(nstars, nlam))

        for istar in range (nstars):
            f.write("{0:e}   {1:e}   {2:e}   {3:e}   {4:e}\n".format(
                radii[istar],
                masses[istar],
                pos[istar, 0],
                pos[istar, 1],
                pos[istar, 2]
            ))

        for ilam in range(nlam):
            f.write("{0:e}\n".format(lam[ilam]))

        for istar in range(nstars):
            if (teff[istar] != 0):
                f.write("{0:f}\n".format(-teff[istar]))
            else:
                for i in range(nlam):
                    f.write("{0:e}\n".format(fstar[ilam]))



def wavelength_micron(lam):
    nlam = len(lam)
    f = open("thermal/wavelength_micron.inp","w")
    f.write("{0:d}\n".format(nlam))
    for ilam in range(nlam):
        f.write("{0:f}\n".format(lam[ilam]))
    f.close()

def mcmono_wavelength_micron(lam_mono):
    nlam = len(lam_mono)

    f = open("thermal/mcmono_wavelength_micron.inp","w")

    f.write("{0:d}\n".format(nlam))
    for ilam in range(nlam):
        f.write("{0:f}\n".format(lam_mono[ilam]))

    f.close()

def amr_grid(radpath, 
             grid, 
             max_level,
             nb_cells, 
             l_cm, 
             gridstyle="regular", 
             coordsystem="cartesian", 
             x=None, 
             y=None, 
             z=None):

    radpath = Path(radpath)
    filepath = radpath / "amr_grid.inp"

    with open(filepath, "w") as f:
        # iformat (typically 1 at present)
        f.write("1\n")

        # grid style 
        if gridstyle == "regular":
            f.write("0\n")
        elif gridstyle == "octtree":
            f.write("1\n")
        elif gridstyle == "amr":
            f.write("10\n")

        # coordinates type
        if coordsystem == "cartesian":
            f.write("0\n")
        elif coordsystem == "spherical":
            f.write("100\n")
        elif coordsystem == "cylindrical":
            f.write("200\n")

        # gridinfo
        f.write("0\n")

        if gridstyle == "octtree":
            f.write("1\t1\t1\n")
            f.write("1\t1\t1\n")
            f.write(f"{max_level}\t{nb_cells}\t{len(grid)}\n")
            for _ in range(3):
                f.write(f"{-l_cm/2:e}\t{l_cm/2:e}\n")
            for g in grid:
                f.write(f"{g}\n")

        elif gridstyle == "regular":
            nx = x.size - 1
            ny = y.size - 1
            nz = z.size - 1

            incl_x = int(nx > 1)
            incl_y = int(ny > 1)
            incl_z = int(nz > 1)

            f.write(f"{incl_x}  {incl_y}  {incl_z}\n")
            f.write(f"{nx}  {ny}  {nz}\n")

            for xi in x:
                f.write(f"{xi:12.9e}\n")
            for yi in y:
                f.write(f"{yi:12.9e}\n")
            for zi in z:
                f.write(f"{zi:12.9e}\n")


        elif (gridstyle == "layer-amr"):
            print("Layer-style AMR grids not yet implemented.")

    # Insert extra info for octtree and amr grids here...

    #f.close()

def dustopac(opacity):
    '''
    Desc: write dustopac.inp
    Args: opacity
    '''
    nspecies = len(opacity)

    f = open("thermal/dustopac.inp","w")

    f.write("2\n")
    f.write("{0:d}\n".format(nspecies))
    f.write("==============================================================\n")
    for i in range(nspecies):
        filetype = opacity[i].split("/")[1].split("_")[0]
        species = opacity[i].split("_")[1].split(".")[0]

        if (filetype == "dustkappa"):
            f.write("1\n")
        elif (filetype == "dustkapscatmat"):
            f.write("10\n")
        elif (filetype == "dustopac"):
            f.write("-1\n")

        f.write("0\n")
        f.write("{0:s}\n".format(species))
        f.write("----------------------------------------------------------\n")

    f.close()


def dust_density(radpath, 
                 density, 
                 nb_cells,
                 nb_sizes,
                 gridstyle="regular"):
    '''
    Desc: write dust_density.inp
    Args: density
    '''
    radpath = Path(radpath)
    filepath = radpath / "dust_density.inp"

    if (gridstyle == "regular"):
        nx, ny, nz = density[0][0].shape
        ncells = nx*ny*nz
        with open(filepath, "w") as f:
            f.write("1\n")
            f.write("{0:d}\n".format(nb_cells))
            f.write("{0:d}\n".format(nb_sizes))

            for ispec in range(len(density[1])): #loop over the dust species within the given structure.
                if (gridstyle == "regular"):
                    for iz in range(nz):
                        for iy in range(ny):
                            for ix in range(nx):
                                f.write("{0:e}\n".format(density[ispec, ix,iy,iz]))

    elif (gridstyle == "octtree"):

        with open(filepath, "w") as f:
            f.write("1\n")
            f.write("{0:d}\n".format(nb_cells))
            f.write("{0:d}\n".format(density.shape[1]))

            for i in range(density.shape[1]):
                for j in range(density.shape[0]):
                    f.write("\n%e" % density[j, i])



def external_rad(isrf):
    '''
    Desc: write external_source.inp
    Args: Interstellar radiation field
    '''
    nlam = len(isrf[0])
    factor = 1e0 #artificially multiply by a factor just to see the impact of different ISRF intensities. Will be removed in future updates
    f = open("thermal/external_source.inp","w")

    f.write("2\n")
    f.write("{0:d}\n".format(nlam))
    for ilam in range(nlam):
        f.write("{0:f}\n".format(isrf[0,ilam]*1e-3))
    for ilam in range(nlam):
        f.write("{0:e}\n".format(isrf[1,ilam]*factor))

    f.close()

def accretion_heating(x, y, z, accretionheating, gridstyle="regular"):
    '''
    Desc: Write viscous accretion heating in the file heatsource.inp 
    Args: Viscous accretion heating structure
    '''
    nx = x.size-1
    ny = y.size-1
    nz = z.size-1
    f = open("thermal/heatsource.inp","w")
    f.write("1\n")
    f.write("{0:d}\n".format(nx*ny*nz))

    if (gridstyle == "regular"):
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    f.write("{0:e}\n".format(accretionheating[ix,iy,iz]))

    f.close()


def numberdens_mol(numberdens, species='CO', gridstyle="regular"):
    '''
    Desc: write numberdens_XXX.inp, where XXX is a chemical species
    Args: 
    '''
    if (gridstyle == "regular"):
        nx, ny = numberdens.shape
        ncells = nx*ny

    print('writing numberdens_{}.inp...'.format(species))

    f = open("thermal/numberdens_{}.inp".format(species),"w")
    f.write("1\n")
    f.write("{0:d}\n".format(ncells))
    if (gridstyle == "regular"):
            for iy in range(ny):
                for ix in range(nx):
                    f.write("{0:e}\n".format(numberdens[ix,iy]))
    f.close()

def lines(species='CO', format='leiden'):
    '''
    Desc: write lines.inp
    Args: species, format
    '''
    f = open("thermal/lines.inp","w")
    f.write("2\n")
    f.write("1\n")
    f.write("{} {} 0 0 0".format(species,format))
    f.close()