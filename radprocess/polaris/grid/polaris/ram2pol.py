# Maire-Anne Carpine
# Ilseung Han
#   29.05.2025
#   04.07.2025

import numpy as np
import math
import struct
import sys
import pymses

from pymses.filters import CellsToPoints
from pymses.utils import constants as C

# constants
grid_id = 20       # grid ID (20 = octree)  
gamma   = 1.4      # Ramses default
mH      = 1.66e-24 # Hydrogen mass in g
kB      = 1.38e-16 # Boltzmann constant in cgs

# X = 0.76 # H mass fraction
# Y = 0.24 # He mass fraction

CLR_LINE =   "                                                      \r"

max_level    = 0
cell_counter = 0
nr_of_cells  = 0

class cell_oct:
    """written by Stefan Reissl"""
    def __init__(self, _x_min, _y_min, _z_min, _length, _level):
        self.x_min = _x_min
        self.y_min = _y_min
        self.z_min = _z_min
        
        self.length = _length
        self.level  = _level
    
        self.isleaf   = 0
        self.data     = []
        self.branches = []  
        
class OcTree:
    """written by Stefan Reissl"""
    def __init__(self, _x_min, _y_min, _z_min, _length):
        self.root = cell_oct(_x_min, _y_min, _z_min, _length, 0)

    def initCellBoundaries(self, cell,_level):
        x_min = cell.x_min
        y_min = cell.y_min
        z_min = cell.z_min
        l     = 0.5 * cell.length

        level = _level

        cell.isleaf      = 0
        cell.data        = []
        cell.branches    = [None, None, None, None, None, None, None, None]
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
            cell_pos.data   = cell.data
            cell_pos.isleaf = 1
                        
        else:    
            if len(cell_pos.branches) == 0:
                self.initCellBoundaries(cell_pos,_level + 1)

            x_mid = cell_pos.x_min + 0.5 * cell_pos.length
            y_mid = cell_pos.y_min + 0.5 * cell_pos.length
            z_mid = cell_pos.z_min + 0.5 * cell_pos.length
            
            new_cell_pos = cell_pos

            if(z_pos < z_mid): #z 0 1 2 3

                if(y_pos < y_mid): #y 0 1

                    if(x_pos < x_mid): #x 0
                        new_cell_pos = cell_pos.branches[0]
                    else: #x 1
                        new_cell_pos = cell_pos.branches[1]

                else: #y 2 3

                    if(x_pos < x_mid): #x 2
                        new_cell_pos = cell_pos.branches[2]
                    else: #x 3
                        new_cell_pos = cell_pos.branches[3]


            else: #z 4 5 6 7

                if(y_pos < y_mid): #y 4 5

                    if(x_pos < x_mid): #x 4
                        new_cell_pos = cell_pos.branches[4]
                    else: #x 5
                        new_cell_pos = cell_pos.branches[5]

                else: #y 6 7

                    if(x_pos < x_mid): #x 6
                        new_cell_pos = cell_pos.branches[6]
                    else: #x 7
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
                
    # def writeOcTree_radmc(self, cell, grid, density):
    #     global cell_counter
    #     global nr_of_cells

    #     if cell.isleaf == 1:    
    #         data_len = len(cell.data)
            
    #         if cell_counter % 10000 == 0:
    #             sys.stdout.write('-> Writing octree grid file : ' + str(100.0 * cell_counter / nr_of_cells) + ' %     \r')
    #             sys.stdout.flush()
                
    #         cell_counter += 1
                
    #         density.append(cell.data[0])
    #         grid.append(0)

    #     else:
    #         grid.append(1)
            
    #         for i in range(8):
    #             self.writeOcTree_radmc(cell.branches[i], grid, density)
        
# def loadRamsesData(filename):
    
#     #Filename will be the format my/folder/here/output_?????
#     #where ????? is from 00001 to 99999
#     print(filename)
    
#     # Split the filename into folder and number (pymses needs this)
#     outstr = "output_"
#     outloc = filename.rfind(outstr)
#     folder = filename[:outloc]
#     numloc = outloc+len(outstr)
    
#     # Note: 5 characters in output number
#     num = int(filename[numloc:numloc+5])
#     # Create the pymses RamsesOutput object 
#     print("\n\nfolder: ", folder , num, "\n\n")
#     snap = pymses.RamsesOutput(folder,num)
    
#     # Create a flat structure with the snapshot's cell data in it
#     print(snap.info)
#     amr = snap.amr_source(["rho","P","vel","B-right","B-left"])
#     print(amr)
#     #amr = snap.amr_source(["rho","P","vel"])
#     cell_source = CellsToPoints(amr)
#     cells = cell_source.flatten()
#     print("The B fields", cells["B-right"],cells["B-left"],"end")

#     # Now make the output dictionary
#     # Spatial information
#     output = {}
    
#     # Cell lengths
#     unit_l = snap.info["unit_length"].express(C.m)
#     output["dx"] = cells.get_sizes()*unit_l
    
#     # max. number of cells
#     numcells = len(output["dx"])
    
#     # Original cell positions (from 0 to 1) converted into uint length)
#     output["x"] = cells.points[:,0]*unit_l
#     output["y"] = cells.points[:,1]*unit_l
#     output["z"] = cells.points[:,2]*unit_l

#     # level of each cell
#     output["level"]=np.log2(unit_l/output["dx"])
    
    
#     # Density in g/cm^3
#     output["dens"] = cells["rho"]*1e6#*snap.info["unit_density"].express(C.g_cc)
#     print ("VAL: DENSITY", output["dens"])
    
#     # Gas temperature in K
#     unit = snap.info["unit_temperature"].express(C.K)
#     X = 0.76 # Hydrogen mass fraction
#     output["Tgas"] = cells["P"]/cells["rho"]*unit/X
    
#     # Velocity in cm/s
#     output["vel"] = cells["vel"]
    
#     # B-field in G
#     unit_b = snap.info['unit_mag'].express(C.Gauss) #(snap.info["unit_pressure"]).express(C.g/C.s/C.s/C.cm)
#     #unit_b = unit_b * 8 * np.pi	
    
#     output["B"] = 0.5*(cells["B-right"]+cells["B-left"])*unit_b
    
#     return output, numcells, unit_l

def loadRamsesData(filename):
    """written by Valeska Valdivia"""
    #Filename will be the format my/folder/here/output_?????
    #where ????? is from 00001 to 99999
    print(filename)
    
    # Split the filename into folder and number (pymses needs this)
    outstr = "output_"
    outloc = filename.rfind(outstr)
    folder = filename[:outloc]
    numloc = outloc+len(outstr)
    
    # Note: 5 characters in output number
    num = int(filename[numloc:numloc+5])
    # Create the pymses RamsesOutput object 
    print("\n\nfolder: ", folder , num, "\n\n")
    snap = pymses.RamsesOutput(folder,num)

    n_dust = snap.info["ndust"] # multi
    
    # Create a flat structure with the snapshot's cell data in it
    print(snap.info)
    dust_names = [f"DR_{i + 1}" for i in range(1, n_dust)] # multi
    amr = snap.amr_source(["rho","P","vel","Br","Bl"]+dust_names)
    print(amr)
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    print("The B fields", cells["Br"], cells["Bl"], "end")
    
    # Now make the output dictionary
    # Spatial information
    output = {}
    
    # Cell lengths
    unit_l = snap.info["unit_length"].express(C.m)
    output["dx"] = cells.get_sizes()*unit_l

    # max. number of cells
    numcells = len(output["dx"])
    
    # Original cell positions (from 0 to 1) converted into uint length)
    output["x"] = cells.points[:,0]*unit_l
    output["y"] = cells.points[:,1]*unit_l
    output["z"] = cells.points[:,2]*unit_l

    print("\nx min max", output["x"].min(), output["x"].max())

    # level of each cell
    output["level"]=np.log2(unit_l/output["dx"])
    
    # Density in cm^-3
    # output["dens"]  = cells["rho"]*1e6
        # this gives it in cm-3
        # *snap.info["unit_density"].express(C.g_cc) this would give the mass density in g/cm-3
    output['denstot'] = cells['rho'] * snap.info['unit_density'].express(C.g_cc)
    output['denstot'] = output['denstot'] * 1e+03
    print("dens min, max", output["denst"].min(), output["denst"].max())

    # Gas temperature in K
    print("*************************************************")
    #============================================
    # Everything in cgs for the moment:
    if 'mu_gas' in snap.info:
        mu=snap.info['mu_gas'] #mass per particles
    else:
        mu=1.4
    print("WORKING WITH mu_gas = ", mu)
    G = 6.7e-8
    kbol  =  1.38062e-16   # erg/degre
    pc=3.08e18 #cm
    mp =  mu * 1.660531e-24  #n gramme
    scale_n = 1.
    scale_l = pc
    scale_d = scale_n*mp
    scale_t = 1.0/np.sqrt(G*scale_d)
    scale_v = scale_l / scale_t    
    scale_T2 = mp/kbol * scale_v**2

    #============================================
    unit = snap.info["unit_temperature"].express(C.K)
    print("unit", unit)
    print("scale_T2", scale_T2)
    X = 0.76 # Hydrogen mass fraction

    output["Tgas"] = cells["P"]/cells["rho"]*scale_T2
    output["Tdust"] = cells["P"]/cells["rho"]*scale_T2
    print("Min Max Mean T", output["Tgas"].min(), output["Tgas"].max(), output["Tgas"].mean())
    print("*************************************************")

    # Velocity in cm/s
    output["vel"] = cells["vel"]*snap.info["unit_velocity"].express(C.cm/C.s)
    
    # B-field in G
    unit_b = snap.info['unit_mag'].express(C.Gauss)
    unit_b = unit_b * 8 * np.pi	
    output["B"] = 0.5*(cells["Br"]+cells["Bl"])*unit_b
    
    # output["mass_dens"] = cells["rho"]*mp # attention, here it is in g/cm3 !!

    # # multi
    # # dust mass densities
    # for i in range(n_dust):
    #     output[f"dust_mass_dens_{i+1}"] = cells[f"DR_{i+1}"]*cells["rho"]*1e6*mp*1e-3  #rho fro cm-3 to m-3, mp from g to kg
    # # enddust
    output['epsilon'] = cells['DR_1'] + cells['DR_2'] + cells['DR_3'] + cells['DR_4'] + cells['DR_5'] +\
                        cells['DR_6'] + cells['DR_7'] + cells['DR_8'] + cells['DR_9'] + cells['DR_10']
    output['densdto'] = output['denstot'] * output['epsilon']
    for i in range(n_dust):
        output[f'densd{(i + 1):02d}'] = output['denstot'] * cells[f'DR_{i + 1}']
    output['densgas'] = output['denstot'] - output['densdto']

    return output, numcells, unit_l, n_dust # multi

def convert_ramses2polaris(datapath, num, outpath, starpos=None, htspot=None,size_hole_au= 4 ):
    global nr_of_cells
    global cell_counter

    o_num    = str(num).zfill(5)
    out_name = "output_"+o_num

    input_file  = datapath+out_name
    output_file = outpath+"ramses_grid_"+o_num+".dat"
    
    print ("input:", input_file)
    print ("output:", output_file)
    
    print ("Loading RAMSES data from: \n", input_file, "\n\n")
    
    # read RAMSES data
    data, nr_of_cells, max_length, n_dust = loadRamsesData(input_file) # multi

    #data IDs for the grid header
    data_ids = [4, 5, 6] + [28, 3] + [29] * n_dust

    # transpose center of cube
    x_min = -0.5*max_length
    y_min = -0.5*max_length
    z_min = -0.5*max_length
    
    max_level = max(data["level"])
    min_level = min(data["level"])
    
    print ("\n\n\n")
    print ("Octree parameter:")
    print ("    Level        (min,max)  : ", int(min_level),",", int(max_level))
    print ("    Nr. of cells (data, max): ", nr_of_cells,",", 8**max_level)
    print ("    Length       (min,max)  : ", max_length/(2**max_level),",", max_length, "\n")
    
    # init. octree
    tree = OcTree(x_min, y_min, z_min, max_length)
            
    #fill octree
    for i in range(0,nr_of_cells):
        
        #create single cell
        level = data["level"][i]
        
        c_x = data["x"][i]-0.5*max_length
        c_y = data["y"][i]-0.5*max_length
        c_z = data["z"][i]-0.5*max_length
        
        # dens = data["dens"][i]
        densgas = data['densgas'][i]

        # # dust
        # dust_dens = np.empty(n_dust, dtype=object)          # multi
        # for j in range(n_dust):                             # multi
        #     dust_dens[j] = data[f"dust_mass_dens_{j+1}"][i] # multi
        # # enddust
        denstot = data['denstot'][i]
        epsilon = data['epsilon'][i]
        densdto = data['densdto'][i]
        densd01 = data['densd01'][i]
        densd02 = data['densd02'][i]
        densd03 = data['densd03'][i]
        densd04 = data['densd04'][i]
        densd05 = data['densd05'][i]
        densd06 = data['densd06'][i]
        densd07 = data['densd07'][i]
        densd08 = data['densd08'][i]
        densd09 = data['densd09'][i]
        densd10 = data['densd10'][i]

        Tgas = data["Tgas"][i]
        
        mag_x = data["B"][i][0]
        mag_y = data["B"][i][1]
        mag_z = data["B"][i][2]
        
        vel_x = data["vel"][i][0]
        vel_y = data["vel"][i][1]
        vel_z = data["vel"][i][2]

        # if starpos is not None:
        #     au_in_m = 1.4959787070e11
        #     dstar = np.sqrt((c_x-starpos[0])**2 + (c_y-starpos[1])**2 + (c_z-starpos[2])**2)  
        #     if (dstar <= size_hole_au*au_in_m):
        #        print("pos", c_x, c_y, c_z)
        #        dens = 0
        #        #dust
        #        for j in range(n_dust):
        #            dust_dens[j] = 0
        #        #enddust
        #        mass_dens = 0
        #        print("dens=", 0,"used to be",data["dens"][i])
        # # END ADDED FROM RADMC2POLARIS
        # # My hotspot 5e14,0,5e14; 6e14,0,5e14 ; 5e14,0,6e14; 5e14,0,7e14 m 
        # # htspot = [[2e14,0,2e14],[3e14,0,2e14],[2e14,0,3e14],[2e14,0,4e14]] #m
        # if htspot is not None:
        #     for hotspot_i in htspot:
        #         au_in_m = 1.4959787070e11
        #         dstar = np.sqrt((c_x-hotspot_i[0])**2 + (c_y-hotspot_i[1])**2 + (c_z-hotspot_i[2])**2)  
        #         if (dstar <= 300.*au_in_m):
        #            print("pos", c_x, c_y, c_z)
        #            dens = 1e8*1e6
        #            #dust
        #            for j in range(n_dust):
        #                dust_dens[j] = 1e8*1e6/n_dust
        #            #enddust
        #            mass_dens = 1e8*1e6*2.3*1.66e-24
        #            print("dens=hotspot","used to be",data["dens"][i])
        
        cell = cell_oct(c_x, c_y, c_z, 0, level)
        # cell.data =[dens, 10, Tgas, mag_x, mag_y, mag_z]
            # XXX here replace dens by all the n_dust dustmassdensities
        # # dust
        # cell.data = dust_dens.tolist() + [dens, 10, Tgas, mag_x, mag_y, mag_z]
            # XXX here replace dens by all the n_dust dustmassdensities
        # # enddust
        # cell.data = [densg, Tgas, 10,
        #              densd1, densd2, densd3, densd4, densd5, densd6, densd7, densd8, densd9, densd10,
        #              mag_x, mag_y, mag_z]
        cell.data = [mag_x, mag_y, mag_z,
                     densgas, Tgas,
                     densd01, densd02, densd03, densd04, densd05, densd06, densd07, densd08, densd09, densd10]
                     # 10]
        
        if i % 10000 == 0:
                sys.stdout.write('Constructing octree: ' + str(100.0 * i / nr_of_cells) + ' %    \r')
                sys.stdout.flush()
        
        #insert single cell into octree
        cell_root= tree.root
        tree.insertInTree(cell_root, cell,0)
    
    sys.stdout.write(CLR_LINE)
    print ("Constructing octree:    done   ")
    
    #check octree integrity
    check = tree.checkOcTree(cell_root)
    sys.stdout.write(CLR_LINE)
        
    if check == False:
        print ("ERROR: Octree integrity is inconsistent!   \n\n")
        exit ()
    else:
        print ("Octree structure   :    OK      ")
    
    #write octree file header
    data_len = len(data_ids)
    file = open(output_file, "wb")
        
    file.write(struct.pack("H", grid_id))
    file.write(struct.pack("H", data_len))

    for d_ids in data_ids:
        file.write(struct.pack("H", d_ids))

    file.write(struct.pack("d", max_length))
    
    #write octree
    cell_counter = 0.0
    tree.writeOcTree(file, tree.root)
    sys.stdout.write(CLR_LINE)

    print ("Writing octree     :    done   \n")
    
    print ("Octree successfully created")
