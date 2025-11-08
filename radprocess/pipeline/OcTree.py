import sys

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
                
    def writeOcTree_radmc(self, cell, grid, density, temp):
        global cell_counter
        global nr_of_cells

        if cell.isleaf == 1:    
            data_len = len(cell.data)
            
            if cell_counter % 10000 == 0:
                sys.stdout.write('-> Writing octree grid file : ' + str(100.0 * cell_counter / nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            cell_counter += 1
                
            density.append(cell.data[0])

            temp.append(cell.data[3])

            grid.append(0)

        else:
            grid.append(1)
            
            for i in range(8):
                self.writeOcTree_radmc(cell.branches[i], grid, density, temp)