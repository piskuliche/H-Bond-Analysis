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
        for icomp in [2,6]:
            self.hyd_by_atom[icomp] = {}
            self.hyd_by_mol[icomp] = {}
        return
    
    def add_hydration_data(self):
        for ifile in range(self.fstart, self.fstop+1):
            print(ifile)
            start_idx = (ifile-1)*self.fr_per_file
            for icomp in [2,6]:
                hyd_mol = self._read_molar_data("%d/hyd_molar_%d.dat"%(ifile,icomp))
                self.hyd_by_mol[icomp] = self._add_arr_to_framewise(self.hyd_by_mol[icomp], hyd_mol, start_idx)
        return
    
    def _read_molar_data(self, filename):
        mol, hyd = np.loadtxt(filename, usecols=(0,1), unpack=True)
        hydration = np.split(hyd, self.fr_per_file)
        return hydration
    
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
    parser.add_argument('-fstart', default=1, type=int, help='Starting index of calculation')
    parser.add_argument('-fstop', default=100, type=int, help='Stopping index of calculation')
    parser.add_argument('-tag', default="5_hyd", type=str, help="Tag for the molecule")
    # Parse the command-line arguments
    args = parser.parse_args()

    hyd_data = Hydration_Data(args.fstart, args.fstop)
    hyd_data.add_hydration_data()
    hyd_data.save_data(args.tag)





