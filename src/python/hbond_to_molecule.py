#!/usr/bin/env python
import numpy as np

class hbond_data:
    def __init__(self, inputname):
        self.inputname = inputname
        self.fr_idx = 0
        self.idx_2_mol = {}
        self.idx_2_comp = {}
        self.components = []
        self.read_hbond_input()
        return
    
    def read_hbond_input(self):
        with open(self.inputname, 'rb') as f:
            for i in range(5): f.readline()
            self.num_components, tmp = f.readline().strip().split()
            self.num_components = int(self.num_components)
            f.readline()
            cnt = 1 # 1 indexing bc fortran
            self.num_atoms = 0
            nmoltmp, apermoltmp = [], []
            for i in range(self.num_components):
                label, nmol, apermol, nacc = f.readline().strip().split()
                check_water = int(f.readline().strip())
                nmol = int(nmol); apermol = int(apermol)
                if check_water == 1:
                    self.nwat = nmol
                nmoltmp.append(nmol); apermoltmp.append(apermol)
                self.num_atoms += nmol*apermol
                for j in range(nmol):
                    for k in range(apermol):
                        self.idx_2_mol[cnt] = j
                        self.idx_2_comp[cnt] = i
                        cnt += 1
            for i in range(self.num_components):
                self.components.append(hbond_component(nmoltmp[i], apermoltmp[i], i))
        return
    
    def read_hbond_data(self, fileno=1):
        mol, h1, h2 = np.loadtxt("%d/all_hydrogen_bonds.dat"%fileno, usecols=(0,1,2), unpack=True)
        h1 = h1.reshape(-1, self.nwat)
        h2 = h2.reshape(-1, self.nwat)
        

        for frame_h1, frame_h2 in zip(h1,h2):
            self.fr_idx += 1
            for component in self.components:
                component.add_frame(self.fr_idx)
            for ev1 in frame_h1[frame_h1!=0]:
                elem = int(ev1)
                component = self.components[self.idx_2_comp[elem]]
                component.add_hbonds(self.fr_idx, self.idx_2_mol[elem])
            for ev2 in frame_h2[frame_h2!=0]:
                elem = int(ev2)
                component = self.components[self.idx_2_comp[elem]]
                component.add_hbonds(self.fr_idx, self.idx_2_mol[elem])
        return

    def generate_hbonds(self, fstart, fstop):
        for i in range(fstart, fstop+1):
            print("File number ", i)
            self.read_hbond_data(i)
    
    def save_map(self, tag):
        import pickle
        savedict = {"idx_2_mol":self.idx_2_mol, "idx_2_comp":self.idx_2_comp}
        pickle.dump(savedict,open("%s_map.pckl"%tag,'wb'))
        return

            
        
class hbond_component:
    def __init__(self, nmol, apermol, cid):
        self.nmol = nmol
        self.apermol = apermol
        self.cid = cid
        self.hbonds = {}

    def add_frame(self, fr_idx):
        self.hbonds[fr_idx] = np.zeros(self.nmol)

    def add_hbonds(self, fr_idx, molec):
        self.hbonds[fr_idx][molec] += 1

    def display(self, fr_idx):
        num_hbonds = np.sum(self.hbonds[fr_idx])
        print("There are %d hbonds in frame %d" % (num_hbonds, fr_idx))
    
    def save_data(self, tag):
        import pickle
        pickle.dump(self.hbonds, open("%s_%d.pckl"%(tag, self.cid),'wb'))


if __name__ == "__main__":
    import argparse
    
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='Choose the options for the analysis setup')

    # Add an argument to the parser
    parser.add_argument('-fstart', default=1, type=int, help='Starting index of calculation')
    parser.add_argument('-fstop', default=100, type=int, help='Stopping index of calculation')
    parser.add_argument('-tag', default="5_hb_", type=str, help="Tag for the molecule")
    parser.add_argument('-map', default = 0, type=int, help="Only generate the map?")
    # Parse the command-line arguments
    args = parser.parse_args()

    data = hbond_data("1/hbonding.in")

    data.save_map(args.tag)

    if args.map == 0:
        data.generate_hbonds(args.fstart,args.fstop)
        for component in data.components:
            component.save_data(args.tag)