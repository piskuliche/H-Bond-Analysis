#!/usr/bin/env python
import numpy as np

class hbond_data:
    def __init__(self, inputname):
        self.inputname = inputname
        self.idx_2_mol = {}
        self.idx_2_comp = {}
        self.components = []
        self.read_hbond_input()
        return
    
    def read_hbond_input(self):
        with open(self.inputname, 'rb') as f:
            for i in range(5): f.readline()
            self.num_components, tmp = f.readline().strip().split()
            f.readline()
            cnt = 1 # 1 indexing bc fortran
            self.num_atoms = 0
            for i in range(self.num_components):
                label, nmol, apermol, nacc = f.readline().strip().split()
                self.components.append(hbond_component(nmol, apermol))
                check_water = int(f.readline().strip())
                if check_water == 1:
                    self.nwat = nmol
                nmol = int(nmol); apermol = int(apermol)
                self.num_atoms += nmol*apermol
                for j in range(nmol):
                    for k in range(apermol):
                        self.idx_2_mol[cnt] = j
                        self.idx_2_comp[cnt] = i
                        cnt += 1
    
    def read_hbond_data(self, fileno=1):
        mol, h1, h2 = np.loadtxt("%d/all_hydrogen_bonds.dat"%fileno, usecols=(0,1,2), unpack=True)
        h1 = h1.reshape(-1, self.nwat)
        h2 = h2.reshape(-1, self.nwat)
        for frame in h1:
            for elem in frame:
                self.components[self.idx_2_comp[elem]].add_hbonds()
            


        
class hbond_component:
    def __init__(self, nmol, apermol, nwat):
        self.nmol = nmol
        self.apermol = apermol
        self.hbonds = {}

    def add_hbonds(self, frame_index):




                        




if __name__ == "__main__":
    run = 1