#!/usr/bin/env python
import numpy as np

class Hydration_Data:
    def __init__(self, fstart=1, fstop=100, num_components=6, fr_per_file=501):
        self.fstart = fstart
        self.fstop = fstop
        self.num_components = num_components
        self.fr_per_file = fr_per_file
        self.hyd_by_atom = {}
        self.hyd_by_mol = {}
        self.hyd_atoms = {}
        self.occ_cmap = {}
        for icomp in [2,6]:
            self.hyd_by_atom[icomp] = {}
            self.hyd_by_mol[icomp] = {}
            self.hyd_atoms[icomp] = {}
            self.occ_cmap[icomp] = {}
        return
    
    def add_hydration_data(self):
        self.add_leaf_data("leaflet_wt10.pckl")
        for ifile in range(self.fstart, self.fstop+1):
            print(ifile)
            start_idx = (ifile-1)*self.fr_per_file
            for icomp in [2,6]:
                hm, hyd_mol = self._read_data("%d/hyd_molar_%d.dat"%(ifile,icomp))
                hyd_atoms, hyd_atom = self._read_data("%d/hyd_atomic_%d.dat"%(ifile,icomp))
                self.hyd_by_mol[icomp] = self._add_arr_to_framewise(self.hyd_by_mol[icomp], hyd_mol, start_idx)
                self.hyd_by_atom[icomp] = self._add_arr_to_framewise(self.hyd_by_atom[icomp], hyd_atom, start_idx)
                self.hyd_atoms[icomp] = self._add_arr_to_framewise(self.hyd_atoms[icomp], hyd_atoms, start_idx)
        return

    def translate_to_mol(self, natoms={2:52, 6:26}):
        all_counts = {}
        for i in [2,6]:
            self.occ_cmap[i] = {0:np.zeros(natoms[i]), 1:np.zeros(natoms[i]), 2:np.zeros(natoms[i]), 3:np.zeros(natoms[i])}
            all_counts[i] = {0:0,1:0,2:0,3:0}
        for frame in self.hyd_by_atom[2]:
            molvals = []
            comps = []
            for atom in self.atomindex[frame]:
                molvals.append(self.map['idx_2_mol'][atom])
                comps.append(self.map['idx_2_comp'][atom]+1)
            for i in range(len(molvals)):
                for atom in range(natoms[comps[i]]):
                    hcount = self.hyd_by_atom[comps[i]][frame][atom]
                    self.occ_cmap[comps[i]][self.occupancy[frame][i]][atom] += hcount
                all_counts[comps[i]][self.occupancy[frame][i]] += 1
        """
        for icomp in [2,6]:
            for occ in self.occ_cmap[icomp]:
                for atom in range(len(self.occ_cmap[icomp][occ])):
                    if all_counts[icomp][occ] > 0:
                        self.occ_cmap[icomp][occ][atom] /= all_counts[icomp][occ]

        """
        for icomp in [2,6]:
            norm = []
            for i in range(4):
                if  np.max(self.occ_cmap[icomp][i]) > 0:
                    norm.append(self.occ_cmap[icomp][i]/np.max(self.occ_cmap[icomp][i]))
                else:
                    norm.append(np.zeros_like(self.occ_cmap[icomp][i]))
            ac = [all_counts[icomp][0], all_counts[icomp][1], all_counts[icomp][2], all_counts[icomp][3]]
            np.savetxt("comp_%d_map.dat"%icomp,np.c_[norm[0], norm[1], norm[2], norm[3]])
            np.savetxt("occ_%d_values.dat"%icomp, np.c_[self.occ_cmap[icomp][0]/ac[0],self.occ_cmap[icomp][1]/ac[1],self.occ_cmap[icomp][2]/ac[2],self.occ_cmap[icomp][3]/ac[3]])

    def add_leaf_data(self,filename):
        import pickle
        self.map, self.atomindex, self.occupancy = pickle.load(open(filename,'rb'))
        return

    
    def _read_data(self, filename):
        mol, hyd = np.loadtxt(filename, usecols=(0,1), unpack=True)
        hydration = np.split(hyd, self.fr_per_file)
        mols = np.split(mol, self.fr_per_file)
        return mols, hydration
    
    def _add_arr_to_framewise(self, framewise, array, start_idx):
        fr_idx = start_idx
        for frame in array:
            framewise[fr_idx] = frame
            fr_idx += 1
        return framewise

    def save_data(self, tag):
        import pickle
        pickle.dump(self.hyd_by_mol, open("%s_mol.pckl"%tag,'wb'))
        return

if __name__ == "__main__":
    import argparse
    
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='Choose the options for the analysis setup')

    # Add an argument to the parser
    parser.add_argument('-fstart', default=50, type=int, help='Starting index of calculation')
    parser.add_argument('-fstop', default=100, type=int, help='Stopping index of calculation')
    parser.add_argument('-tag', default="5_hyd", type=str, help="Tag for the molecule")
    # Parse the command-line arguments
    args = parser.parse_args()
    hyd_data = Hydration_Data(args.fstart, args.fstop)
    hyd_data.add_hydration_data()
    #hyd_data.save_data(args.tag)
    hyd_data.translate_to_mol()






