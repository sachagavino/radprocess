import sys
import struct

class CellOct:
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
        self.root = CellOct(_x_min, _y_min, _z_min, _length, 0)
        self.cell_counter = 0
        self.nr_of_cells = 0

    def reset_counter(self):
        """Reset the cell counter to zero."""
        self.cell_counter = 0



    def initCellBoundaries(self, cell,_level):
        x_min = cell.x_min
        y_min = cell.y_min
        z_min = cell.z_min
        l = 0.5 * cell.length

        level = _level

        cell.isleaf = 0
        cell.data = []
        cell.branches = [None, None, None, None, None, None, None, None]
        cell.branches[0] = CellOct(x_min, y_min, z_min, l, level)
        cell.branches[1] = CellOct(x_min + l, y_min, z_min, l, level)
        cell.branches[2] = CellOct(x_min, y_min + l, z_min, l, level)
        cell.branches[3] = CellOct(x_min + l, y_min + l, z_min, l, level)

        cell.branches[4] = CellOct(x_min, y_min, z_min + l, l, level)
        cell.branches[5] = CellOct(x_min + l, y_min, z_min + l, l, level)
        cell.branches[6] = CellOct(x_min, y_min + l, z_min + l, l, level)
        cell.branches[7] = CellOct(x_min + l, y_min + l, z_min + l, l, level)     
        
    def insertInTree(self, cell_pos, cell, _level):    
        x_pos = cell.x_min
        y_pos = cell.y_min
        z_pos = cell.z_min

        if cell_pos.level == cell.level:
            cell_pos.data=cell.data
            cell_pos.isleaf=1
                        
        else:    
            if len(cell_pos.branches)==0:
                self.initCellBoundaries(cell_pos,_level+1)

            x_mid = cell_pos.x_min+0.5*cell_pos.length
            y_mid = cell_pos.y_min+0.5*cell_pos.length
            z_mid = cell_pos.z_min+0.5*cell_pos.length
            
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
        """
        Write the octree to a binary file in POLARIS format.
        
        Handles the case where branch nodes have children that were
        never populated (e.g., in subbox extraction where not all 8
        octants contain cells). Unpopulated children are written as
        empty leaf cells.
        """
        # Handle None cells (should not happen, but defensive)
        if cell is None:
            return

        # If this is a branch node but has no children, treat it as a leaf
        # This happens in subbox extraction when an octant has no cells
        if cell.isleaf == 0 and len(cell.branches) == 0:
            cell.isleaf = 1
            if len(cell.data) == 0:
                # Write as empty leaf — need to know data length
                # Use _n_data if available, otherwise skip
                n_data = getattr(self, '_n_data', 0)
                cell.data = [0.0] * n_data

        file.write(struct.pack("H", cell.isleaf))
        file.write(struct.pack("H", cell.level))   

        if cell.isleaf == 1:    
            data_len = len(cell.data)
            
            if self.cell_counter % 10000 == 0:
                sys.stdout.write('-> Writing octree grid file : ' + str(100.0 * self.cell_counter / self.nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            self.cell_counter += 1 
         
            for i in range(0, data_len):
                file.write(struct.pack("f", cell.data[i]))
        else:
            for i in range(8):
                self.writeOcTree(file, cell.branches[i])
                
                
    def checkOcTree(self, cell):

        if cell is None:
            return True

        if cell.isleaf == 1:    
            length = len(cell.data)
            
            if length == 0:
                return False
            
            
            if self.cell_counter % 10000 == 0:
                sys.stdout.write('-> Checking octree integrity : ' + str(100.0 * self.cell_counter / self.nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            self.cell_counter += 1    
            
        else:
            length = len(cell.branches)
            
            if length == 0:
                # Empty branch in subbox — this is OK, will be written as empty leaf
                return True
            
            for i in range(8):
                if not self.checkOcTree(cell.branches[i]):
                    return False
                
        return True
                
    def writeOcTree_radmc(self, cell, grid, density):
        if cell is None:
            # Should not happen, but safety
            grid.append(0)
            density.append([0.0] * self._n_species)
            self.cell_counter += 1
            return

        if cell.isleaf == 1:    
            if self.cell_counter % 10000 == 0:
                sys.stdout.write('-> Writing octree grid file : ' + str(100.0 * self.cell_counter / self.nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            self.cell_counter += 1
            density.append(cell.data)
            grid.append(0)

        else:
            if len(cell.branches) == 0:
                # Empty intermediate node — treat as leaf with zero density
                if self.cell_counter % 10000 == 0:
                    sys.stdout.write('-> Writing octree grid file : ' + str(100.0 * self.cell_counter / self.nr_of_cells) + ' %     \r')
                    sys.stdout.flush()
                self.cell_counter += 1
                grid.append(0)
                density.append([0.0] * self._n_species)
            else:
                grid.append(1)
                for i in range(8):
                    self.writeOcTree_radmc(cell.branches[i], grid, density)