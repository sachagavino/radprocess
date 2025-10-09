# Marie-Anne Carpine
# Ilseung Han
#   27.07.2025
#   31.07.2025
#   13.08.2025

import os
import sys
import math
import pymses
import struct
# Import modules
import numpy as np
import datetime as dt

from pymses.utils import constants as C
from pymses.filters import CellsToPoints

# constants
grid_id = 20         # grid ID (20 = octree)  
gamma   = 1.4        # Ramses default
# mH    = 1.6726e-24 # Hydrogen mass in g # 25.03.2025 (commented out)
mH      = 1.66e-24   # Hydrogen mass in g # 25.03.2025 (added)
    # the same as converter_ramses.py
kB      = 1.38e-16   # Boltzmann constant in cgs

CLR_LINE =   "                                                      \r"

max_level    = 0
cell_counter = 0
nr_of_cells  = 0

class cell_oct:
    """written by Stefan Reiss"""
    def __init__(self, _x_min, _y_min, _z_min, _length, _level):
        self.x_min = _x_min
        self.y_min = _y_min
        self.z_min = _z_min
        
        self.length = _length
        self.level = _level
    
        self.isleaf = 0
        self.data = []
        self.branches = []  
                  
class OcTree:
    """written by Stefan Reiss"""
    def __init__(self, _x_min, _y_min, _z_min, _length):
        self.root = cell_oct(_x_min, _y_min, _z_min, _length, 0)

    def initCellBoundaries(self, cell,_level):
        x_min = cell.x_min
        y_min = cell.y_min
        z_min = cell.z_min
        l = 0.5 * cell.length

        level = _level

        cell.isleaf = 0
        cell.data = []
        cell.branches = [None, None, None, None, None, None, None, None]
        cell.branches[0] = cell_oct(x_min, y_min, z_min, l, level)
        cell.branches[1] = cell_oct(x_min + l, y_min, z_min, l, level)
        cell.branches[2] = cell_oct(x_min, y_min + l, z_min, l, level)
        cell.branches[3] = cell_oct(x_min + l, y_min + l, z_min, l, level)

        cell.branches[4] = cell_oct(x_min, y_min, z_min + l, l, level)
        cell.branches[5] = cell_oct(x_min + l, y_min, z_min + l, l, level)
        cell.branches[6] = cell_oct(x_min, y_min + l, z_min + l, l, level)
        cell.branches[7] = cell_oct(x_min + l, y_min + l, z_min + l, l, level)     
        
    def insertInTree(self, cell_pos, cell, _level):    
        x_pos = cell.x_min
        y_pos = cell.y_min
        z_pos = cell.z_min

        if cell_pos.level == cell.level:
            cell_pos.data = cell.data
            cell_pos.isleaf = 1
                        
        else:    
            if len(cell_pos.branches) == 0:
                self.initCellBoundaries(cell_pos, _level + 1)

            x_mid = cell_pos.x_min + 0.5 * cell_pos.length
            y_mid = cell_pos.y_min + 0.5 * cell_pos.length
            z_mid = cell_pos.z_min + 0.5 * cell_pos.length
            
            new_cell_pos = cell_pos

            if(z_pos < z_mid): # z 0 1 2 3

                if(y_pos < y_mid): # y 0 1

                    if(x_pos < x_mid): # x 0
                        new_cell_pos = cell_pos.branches[0]
                    else: # x 1
                        new_cell_pos = cell_pos.branches[1]

                else: # y 2 3

                    if(x_pos < x_mid): # x 2
                        new_cell_pos = cell_pos.branches[2]
                    else: # x 3
                        new_cell_pos = cell_pos.branches[3]

            else: # z 4 5 6 7

                if(y_pos < y_mid): # y 4 5

                    if(x_pos < x_mid): # x 4
                        new_cell_pos = cell_pos.branches[4]
                    else: # x 5
                        new_cell_pos = cell_pos.branches[5]

                else: # y 6 7

                    if(x_pos < x_mid): # x 6
                        new_cell_pos = cell_pos.branches[6]
                    else: # x 7
                        new_cell_pos = cell_pos.branches[7]

            self.insertInTree(new_cell_pos, cell, _level+1)


    def writeOcTree(self, file, cell):
        global cell_counter
        global nr_of_cells

        file.write(struct.pack("H", cell.isleaf))
        file.write(struct.pack("H", cell.level))   

        if cell.isleaf == 1:    
            data_len = len(cell.data)
            
            if cell_counter % 10000 == 0:
                sys.stdout.write('-> Writing octree grid file : ' + str(100.0 * cell_counter / nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            cell_counter += 1 
         
            for i in range(0, data_len):
                file.write(struct.pack("f", cell.data[i]))
        else:
            for i in range(8):
                self.writeOcTree(file, cell.branches[i])
                
                
    def checkOcTree(self, cell):
        global cell_counter
        global nr_of_cells

        if cell.isleaf == 1:    
            length = len(cell.data)
            
            if length == 0:
                return False
            
            if cell_counter % 10000 == 0:
                sys.stdout.write('-> Checking octree integrity : ' + str(100.0 * cell_counter / nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            cell_counter += 1    
            
        else:
            length = len(cell.branches)
            
            if length == 0:
                return False
            
            for i in range(8):
                self.checkOcTree(cell.branches[i])                
                
        return True

    # RADMC-3D                
    def writeOcTree_radmc(self, cell, grid, density):
        global cell_counter
        global nr_of_cells

        if cell.isleaf == 1:    
            data_len = len(cell.data)
            
            if cell_counter % 10000 == 0:
                sys.stdout.write('-> Writing octree grid file : ' + str(100.0 * cell_counter / nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            cell_counter += 1
                
            # density.append(cell.data[0])
            # density.append(cell.data[1])
            density.append(cell.data)
            grid.append(0)

        else:
            grid.append(1)
            
            for i in range(8):
                self.writeOcTree_radmc(cell.branches[i], grid, density)
        
def loadRamsesData(filename):
    """written by Valeska Valdivia"""
    # Filename will be the format my/folder/here/output_?????
    # where ????? is from 00001 to 99999
    print(filename)
    
    # Split the filename into folder and number (pymses needs this)
    outstr = "output_"
    outloc = filename.rfind(outstr)
    folder = filename[:outloc]
    numloc = outloc + len(outstr)
    
    # Note: 5 characters in output number
    num = int(filename[numloc:numloc+5])
    # Create the pymses RamsesOutput object 
    print("\n\nfolder: ", folder , num, "\n\n")
    snap = pymses.RamsesOutput(folder,num)

    n_dust = snap.info["ndust"] # multi

    # Create a flat structure with the snapshot's cell data in it
    print(snap.info)
    snap.amr_fields()
    dust_names = [f"DR_{i + 1}" for i in range(n_dust)] # multi
    amr = snap.amr_source(["rho", "P", "vel", "Br", "Bl"] + dust_names)
    print(amr)
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    print("The B fields", cells["Br"], cells["Bl"], "end")
    
    # Now make the output dictionary
    # Spatial information
    output = {}
    
    # Cell lengths
    unit_l = snap.info["unit_length"].express(C.m)
    output["dx"] = cells.get_sizes() * unit_l

    # max. number of cells
    numcells = len(output["dx"])
    
    # Original cell positions (from 0 to 1) converted into uint length)
    output["x"] = cells.points[:,0] * unit_l
    output["y"] = cells.points[:,1] * unit_l
    output["z"] = cells.points[:,2] * unit_l

    print("\nx min max", output["x"].min(), output["x"].max())

    # level of each cell
    output["level"] = np.log2(unit_l / output["dx"])
    
    # # number density in cgs
    # output["dens"]  = cells["rho"]  # this gives it in cm-3
        # (25.03.2025) This part is different from 'converter_ramses.py'
        # *snap.info["unit_density"].express(C.g_cc) this would give the mass density in g/cm-3
    # mass density in cgs
    output['denstot'] = cells['rho'] * snap.info['unit_density'].express(C.g_cc)
    print("dens min, max", output["denstot"].min(), output["denstot"].max())

    # Gas temperature in K
    print("*************************************************")
    #============================================
    # Everything in cgs for the moment:
    if 'mu_gas' in snap.info:
        mu = snap.info['mu_gas'] # mass per particles
    else:
        mu = 1.4
    print("WORKING WITH mu_gas = ", mu)
    G        = 6.7e-8
    kbol     = 1.38062e-16          # erg/degre
    pc       = 3.08e18              # cm
    mp       = mu * 1.660531e-24    # n gramme
    scale_n  = 1.
    scale_l  = pc
    scale_d  = scale_n * mp
    scale_t  = 1.0 / np.sqrt(G * scale_d)
    scale_v  = scale_l / scale_t    
    scale_T2 = mp/kbol * scale_v**2

    #============================================
    unit = snap.info["unit_temperature"].express(C.K)
    print("unit", unit)
    print("scale_T2", scale_T2)
    X = 0.76 # Hydrogen mass fraction

    output["Tgas"] = cells["P"] / cells["rho"] * scale_T2
    output["Tdust"] = cells["P"] / cells["rho"] * scale_T2
    print("Min Max Mean T", output["Tgas"].min(), output["Tgas"].max(), output["Tgas"].mean())
    print("*************************************************")

    # Velocity in m/s
    output["vel"] = cells["vel"] * snap.info["unit_velocity"].express(C.m / C.s)
    
    # B-field in G
    unit_b      = snap.info['unit_mag'].express(C.T)
    unit_b      = unit_b * 8 * np.pi
    output["B"] = 0.5 * (cells["Br"] + cells["Bl"]) * unit_b
   
    # mass density in cgs
    # output["mass_dens"] = cells["rho"] * mp

    # # multi
    output['epsilon'] = cells['DR_1'] + cells['DR_2'] + cells['DR_3'] + cells['DR_4'] + cells['DR_5'] +\
                        cells['DR_6'] + cells['DR_7'] + cells['DR_8'] + cells['DR_9'] + cells['DR_10']
    output['densdto'] = output['denstot'] * output['epsilon']
    for i in range(n_dust):
        output[f'densd{(i + 1):02d}'] = output['denstot'] * cells[f'DR_{i + 1}']
    output['densgas'] = output['denstot'] - output['densdto']

    # return output, numcells, unit_l
    return output, numcells, unit_l, n_dust # multi

# def convert_ramses2radmc(datapath, num, dust_to_gas, outpath, starpos = None, htspot = None, size_hole_au = 4):
def convert_ramses2radmc(datapath, num, outpath, starpos = None, htspot = None, size_hole_au = 4):
    """
    Function to call loadRamsesData() to load RAMSES outputs and write 'dust_density.inp' and 'amr_grid.inp' for radmc3d

    INPUTS:
    datapath (string)   : Path to the simulations data (that contains different output folders)
    num (int)           : Simulation output number
    dust_to_gas (float) : Dust to gas mass ratio
    outpath (string)    : Directory to write output files to

    OUTPUT: 'dust_density.inp' and 'amr_grid.inp' files in outpath
    """
    global nr_of_cells
    global cell_counter
   
    print(pymses.__file__)

    print ("====================================================================================================")
    o_num      = str(num).zfill(5)
    out_name   = "output_" + o_num
    input_file = datapath + out_name

    print ("Ramses input:    ", input_file)
    print ("Output for RADMC: amr_grid.inp, dust_density.inp")
    print ("")

    print ("Loading RAMSES data from: \n", input_file)
    #
    # Read RAMSES data
    #
    # data, nr_of_cells, max_length = loadRamsesData(input_file)
    data, nr_of_cells, max_length, n_dust = loadRamsesData(input_file) # multi
    # L_cm = pymses.RamsesOutput(datapath, num).info["unit_length"].express(C.cm) # box length in cm
    L_cm = pymses.RamsesOutput(datapath, int(num)).info["unit_length"].express(C.cm) # box length in cm
    #
    # Transpose center of cube
    #
    x_min = -0.5 * max_length
    y_min = -0.5 * max_length
    z_min = -0.5 * max_length
    
    max_level = max(data["level"])
    min_level = min(data["level"])
    
    print("\n")
    print("Octree parameter:")
    print("    Level        (min,max)  : ", int(min_level), ",", int(max_level))
    print("    Nr. of cells (data, max): ", nr_of_cells, ",", 8**max_level)
    print("    Length       (min,max)  : ", max_length / (2**max_level), ",", max_length, "\n")
    #
    # Init. octree
    #
    tree = OcTree(x_min, y_min, z_min, max_length)
    #
    # Fill octree
    #
    for i in range(0, nr_of_cells):
        #
        # Create single cell
        #
        level = data["level"][i]
        
        c_x = data["x"][i] - 0.5 * max_length
        c_y = data["y"][i] - 0.5 * max_length
        c_z = data["z"][i] - 0.5 * max_length
        
        # mass_dens = data["mass_dens"][i]
        # dens      = data["dens"][i]
        # dust_dens = np.empty(n_dust, dtype = object)         # multi
        # for j in range(n_dust):                              # multi
        #     dust_dens[j] = data[f"dust_mass_dens{j + 1}"][i] # multi
        denstot   = data['denstot'][i]
        epsilon   = data['epsilon'][i]
        densdto   = data['densdto'][i]
        densd01   = data['densd01'][i]
        densd02   = data['densd02'][i]
        densd03   = data['densd03'][i]
        densd04   = data['densd04'][i]
        densd05   = data['densd05'][i]
        densd06   = data['densd06'][i]
        densd07   = data['densd07'][i]
        densd08   = data['densd08'][i]
        densd09   = data['densd09'][i]
        densd10   = data['densd10'][i]
        # Tgas      = data["Tgas"][i]
        # Tdust     = data["Tdust"][i]
        
        # mag_x = data["B"][i][0]
        # mag_y = data["B"][i][1]
        # mag_z = data["B"][i][2]
        
        # vel_x = data["vel"][i][0]
        # vel_y = data["vel"][i][1]
        # vel_z = data["vel"][i][2]
        
        if starpos is not None:
            au_in_m = 1.4959787070e11
            dstar = np.sqrt((c_x - starpos[0])**2 + (c_y - starpos[1])**2 + (c_z - starpos[2])**2)  
            if (dstar <= size_hole_au * au_in_m):
                print("pos", c_x, c_y, c_z)
                # dens      = 0
                # mass_dens = 0
                # print("dens=",  0, "used to be", data["dens"][i])
                densd01 = 0
                densd02 = 0
                densd03 = 0
                densd04 = 0
                densd05 = 0
                densd06 = 0
                densd07 = 0
                densd08 = 0
                densd09 = 0
                densd10 = 0
                print('densd01 = 0 used to be', data['densd01'][i])

        # # My hotspot 5e14,0,5e14; 6e14,0,5e14 ; 5e14,0,6e14; 5e14,0,7e14 m 
        # # htspot = [[2e14,0,2e14],[3e14,0,2e14],[2e14,0,3e14],[2e14,0,4e14]] # m
        # if htspot is not None:
        #     for hotspot_i in htspot:
        #         au_in_m = 1.4959787070e11
        #         dstar = np.sqrt((c_x - hotspot_i[0])**2 + (c_y - hotspot_i[1])**2 + (c_z - hotspot_i[2])**2)  
        #         if (dstar <= 300. * au_in_m):
        #            print("pos", c_x, c_y, c_z)
        #            dens      = 1e8
        #            mass_dens = 1e8 * 2.3 * 1.66e-24
        #            print("dens=hotspot", "used to be", data["dens"][i])
        # # END HSTPOT
        
        cell = cell_oct(c_x, c_y, c_z, 0, level)
        # cell.data = [mass_dens, dens, Tdust, Tgas, mag_x, mag_y, mag_z, vel_x, vel_y, vel_z, 0]
        #     # for single dust distribution
        cell.data = [densd01, densd02, densd03, densd04, densd05,
                     densd06, densd07, densd08, densd09, densd10]
        #
        # Insert single cell into octree
        #
        cell_root = tree.root
        tree.insertInTree(cell_root, cell, 0)

        if i % 10000 == 0:
            sys.stdout.write('Constructing octree: ' + str(100.0 * i / nr_of_cells) + ' %    \r')
            sys.stdout.flush()
    #
    #######
    #print("Extrema")
    #print(c_x_min,c_y_min,c_z_min)
    #print(c_x_max,c_y_max,c_z_max)
    #######
    sys.stdout.write(CLR_LINE)
    print ("Constructing octree:    done   ")
    #
    # Check octree integrity
    #
    print ("Calling tree.checkOcTree(cell_root), nr_of_cells=", nr_of_cells)
    check = tree.checkOcTree(cell_root)
    print ("Tree OK ;)")
    sys.stdout.write(CLR_LINE)
        
    if check == False:
        print ("ERROR: Octree integrity is inconsistent!   \n\n")
        exit ()
    else:
        print ("Octree structure   :    OK      ") 
            
    print("Writing octree     :    done   \n")
    
    print("Octree successfully created\n")
    #
    # Write octree
    #
    cell_counter = 0.0    
    grid         = []   # Vector containing cell values (0/1) in RADMC convention
    density      = []   # Vector containing dust density in RADMC convention
    tree.writeOcTree_radmc(tree.root, grid, density)
    densityarray = np.array(density)
    sys.stdout.write(CLR_LINE)
    #
    # Create output directory if not existed
    #
    os.makedirs(os.path.dirname(outpath), exist_ok = True)
    #
    print("Writing the amr_grid.inp file for RADMC-3D...\n")
    #
    with open(outpath+'amr_grid.inp','w+') as f:
        f.write('1\n')                                                  # iformat <=== Typically 1 at present
        f.write('1\n')                                                  # Grid style (1 = Oct-tree)
        f.write('1\n')                                                  # coordsystem (cartesian if coordsystem < 100)
        f.write('0\n')                                                  # gridinfo (= 0 recommended)
        f.write('1\t1\t1\n')                                            # incl_x,incl_y,incl_z
        f.write('1\t1\t1\n')                                            # nx,ny,nz
        f.write('%d\t%d\t%d\n' % (max_level, nr_of_cells, len(grid)))   # levelmax, nleafsmax, nbranchmax
        f.write('%e\t%e\n' % (-L_cm / 2, L_cm / 2))
        f.write('%e\t%e\n' % (-L_cm / 2, L_cm / 2))
        f.write('%e\t%e\n' % (-L_cm / 2, L_cm / 2))

        for i in range(len(grid)):
            f.write('%d\n' % grid[i])
    #
    print("Writing the dust_density.inp file for RADMC-3D...\n")
    #
    with open(outpath + 'dust_density.inp','w+') as f:
        f.write('1\n')
        f.write('%d\n' % nr_of_cells)
        f.write('%d' % densityarray.shape[1])
       
        for i in range(densityarray.shape[1]):
            for j in range(densityarray.shape[0]):
                f.write('\n%e' % densityarray[j, i])
        # return density, densityarray
