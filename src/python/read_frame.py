#!/usr/bin/env python
import numpy as np
import pickle

""" 
********************************************
CLASS: FRAME_DATA

This class stores the data for a single frame of the calculation.
It assumes that you have stored hydration data, hydrogen bonding data, 
voronoi data, as well as lipid/laurdan conformation data in the appropriate formats (more the other class)
********************************************
"""
class frame_data():

    def __init__(self, nmols, ncomps=[180,20], compno=[2,6]):
        self.components = [2] * ncomps[0] + [6] * ncomps[1]
        self.nmols = nmols
        
        # Voronoi Data
        self.data = {}
        self.data["hyd_mol"] = np.zeros(nmols)
        self.data['areas'] = np.zeros(nmols)
        self.data['occupancy'] = np.zeros(nmols)
        self.data['leaflet'] = np.zeros(nmols)
        self.data['molids'] = np.zeros(nmols)
        self.notes = []

        self.compno = compno
        self.ncomps = ncomps
        return
    
    def add_data(self, dataname, data):
        """ Adds data to the frame.

        Args:
            dataname (str): The name of the data to be added.
            data (float(s)): The data to be added.
        
        Returns:
            None
        """
        if dataname not in self.data:
            self.data[dataname] = np.zeros(self.nmols)
        for i, d in enumerate(data):
            self.data[dataname][i] = d
        return 
    
    def grab_data(self, dataname, occupancy=False):
        """ Returns the data from the frame.

        Args:
            dataname (str): The name of the data to be returned.
            occupancy (bool, optional): Whether or not to return the data by occupancy. Defaults to False.

        Returns:
            array: The data from the frame - split by mol

        """
        if occupancy:
            if "occupancy" not in self.data:
                raise ValueError("Occupancy data not present in frame.")
            # Masks the data by occupancy
            occdata0 = self.data[dataname][(self.data['occupancy'] == 0)]
            occdata1 = self.data[dataname][(self.data['occupancy'] == 1)]
            occdata2 = self.data[dataname][(self.data['occupancy'] == 2)]
            occdata3 = self.data[dataname][(self.data['occupancy'] == 3)]
            return occdata0, occdata1, occdata2, occdata3
        else:
            return self.data[dataname]
        
    
    def query_average(self, dataname, molec=2, occupancy=False):
        """ Returns the average of the data in the frame.

        Args:
            dataname (str): The name of the data to be averaged.
            occupancy (bool, optional): Whether or not to average the data by occupancy. Defaults to False.

        Returns:
            float(s): The average of the data.

        """
        if occupancy:
            data = self.grab_data(dataname, molec=molec, occupancy=True)
            avs = []
            for i in range(4):
                if len(data[i]) > 0:
                    avs.append(np.average(data[i]))
                else:
                    avs.append(0)
                
            return avs
        else:
            return np.average(self.grab_data(dataname, molec=molec))
        
    def histogram_frame(self, dataname, molec=2, bins=50, range=(0,100), 
                        bin_edges = None, occupancy=False):
        """
        """
        data = self.grab_data(dataname, occupancy=occupancy)[molec]
        if bin_edges is None:
            hist, bin_edges = np.histogram(data, bins=bins, range=range)
        else:
            hist, _ = np.histogram(data, bins=bin_edges)
        return hist, bin_edges
    
    def add_note(self, note):
        """ Adds a note to the frame.

        Args:
            note (str): The note to be added.

        Returns:
            None

        """
        self.notes.append(note)
        return


""" 
********************************************
CLASS: CALCULATION_DATA

This class builds and pulls data to fill out frame_data objects in order to make it possible
to decompose the objects by occupancy etc. It also contains some basic functionality for analysis
though that is more meant to be plugged into other modules.

The key dependency is that it requires access to the frame_data class.

********************************************
"""

class calculation_data:

    def __init__(self, dirframes=501, dirstart=1, dirstop=100, ncomps=[180,2], compno=[2,6]):
        """ Initializes the class
        """
        self.dirframes = dirframes
        self.dirstart  = dirstart
        self.dirstop   = dirstop
        self.ncomps    = ncomps
        self.compno    = compno
        self.nmols     = sum(ncomps)
        self.frames    = []
        self.build_data()
        return
    
    def build_data(self):
        """ Builds the data from the directories.

        Args:
            None

        Returns:
            None
        """
        for i in range(self.dirstart, self.dirstop):
            self.new_dir()
            self.pull_dir_hyd(i)
            self.pull_voronoi_data(i)
            self.pull_lipid_data(i)
        return

    def new_dir(self):
        """ Creates a new set of frames in the data

        Args:
            None
        
        Returns:
            None

        """
        for frame in range(self.dirframes):
            self.frames.append(frame_data(self.nmols, self.ncomps, self.compno))
        return
    
    def pull_lipid_data(self, ndir):
        """ Pulls data from the lipid calculations

        Args:

        Returns:

        """
        zdists = pickle.load(open("%d/output_zdists_%d.pckl" % (ndir,ndir), 'rb'))
        rg = pickle.load(open("%d/output_rg_%d.pckl" % (ndir,ndir), 'rb'))
        p2 = pickle.load(open("%d/output_laurp2_%d.pckl" % (ndir,ndir), 'rb'))
        for frame in range(self.dirframes):
            self.frames[frame+(ndir-1)*self.dirframes].add_data("zdists", zdists[frame])
            self.frames[frame+(ndir-1)*self.dirframes].add_data("rg", rg[frame])
            self.frames[frame+(ndir-1)*self.dirframes].add_data("laurtilt", p2[frame])
        return

    def pull_dir_hyd(self, ndir):
        """ Pulls the hyd data from the directory.

        Args:
            ndir (int): The directory number to pull the data from.

        Returns:
            None

        """
        hyd_comp = {}

        for comp in self.compno:
            hyd_comp[comp] = self.hyd_molar_data(ndir, comp)

        for frame in range(self.dirframes):
            hyd_data = np.array([])
            for comp in self.compno:
                hyd_data = np.append(hyd_data, hyd_comp[comp][frame])
            self.frames[frame+(ndir-1)*self.dirframes].add_data("hyd_mol", hyd_data)
        return
    
    def pull_voronoi_data(self, ndir):
        """The Voronoi data is stored as a pickle file which is segregated by frame and leaf.

        voronoi pickles include:

        occupancy
        areas
        resnames
        molids
        atomindex

        Args:
            ndir (int): The directory number to pull the data from.

        """
        try:
            vor_data = pickle.load(open('../voronoi_plots/data/vor_%d.pckl' % ndir, 'rb'))
        except FileNotFoundError:
            print("No voronoi data for directory %d" % ndir)
            for frame in range(self.dirframes):
                frame_index = frame + (ndir-1) * self.dirframes
                mols = np.sum(self.ncomps)
                self.frames[frame_index].add_data("occupancy", np.zeros(mols))
                self.frames[frame_index].add_data("areas", np.zeros(mols))
                self.frames[frame_index].add_data("leaflet", np.zeros(mols))
                self.frames[frame_index].add_data("molids", np.zeros(mols))
                self.frames[frame_index].add_note("voronoi: No voronoi data for directory %d" % ndir)
            return
        for frame in range(self.dirframes):
            # Choose the correct frame
            frame_index = frame + (ndir-1) * self.dirframes

            # Pull the data from the pickle
            print("fr",frame, ndir)
            leaf1 = vor_data[frame_index][0]
            leaf2 = vor_data[frame_index][1]
            
            # Geneate the leaflet array 0 for one leaflet, 1 for the other.
            l1s, l2s = np.zeros(len(leaf1['molids'])), np.ones(len(leaf2['molids']))

            # Concatenate the data into a single array
            molids  = np.concatenate((leaf1['molids'],leaf2['molids']))
            leafs   = np.concatenate((l1s,l2s))
            areas   = np.concatenate((leaf1['areas'],leaf2['areas']))
            occ     = np.concatenate((leaf1['occupancy'],leaf2['occupancy']))

            # Sort the data by molecule id
            sort_indices = np.argsort(molids)

            # Add the data to the frame
            self.frames[frame_index].add_data("occupancy", occ[sort_indices])
            self.frames[frame_index].add_data("areas", areas[sort_indices])
            self.frames[frame_index].add_data("leaflet", leafs[sort_indices])
            self.frames[frame_index].add_data("molids", molids[sort_indices])


        return
    
    def hyd_molar_data(self, ndir, comp):
        """ Pulls the hyd molar data from the directory.

        Args:
            ndir (int): The directory number to pull the data from.
            comp (int): The component number to pull the data from.

        Returns:
            np.array: The hyd molar data for the component.

        """
        hyd = np.loadtxt('%d/hyd_molar_%d.dat' % (ndir, comp) , usecols=(1), unpack=True)
        return np.array_split(hyd, self.dirframes)
    
    def analyze_data(self, dataname, start=None, stop=None, occupancy=False):
        if start == None:
            start = self.dirstart
        if stop == None:
            stop = self.dirstop

    
    def grab_frame_data(self, dataname, start, stop, molec=2, occupancy=False, flatten=False, average=False):
        if flatten == True and average == True:
            raise ValueError("Cannot flatten and average data.")
        if occupancy:
            data1, data2, data3, data4 = [], [], [], []
            for frame in range(start, stop):
                data1.append(self.frames[frame].grab_data(dataname, molec=molec, occupancy=True)[0])
                data2.append(self.frames[frame].grab_data(dataname, molec=molec, occupancy=True)[1])
                data3.append(self.frames[frame].grab_data(dataname, molec=molec, occupancy=True)[2])
                data4.append(self.frames[frame].grab_data(dataname, molec=molec, occupancy=True)[3])
            if flatten:
                # flatten arrays
                data1 = [item for sublist in data1 for item in sublist]
                data2 = [item for sublist in data2 for item in sublist]
                data3 = [item for sublist in data3 for item in sublist]
                data4 = [item for sublist in data4 for item in sublist]
            if average:
                # average arrays
                data1 = np.mean(data1, axis=1)
                data2 = np.mean(data2, axis=1)
                data3 = np.mean(data3, axis=1)
                data4 = np.mean(data4, axis=1)
            return np.array(data1), np.array(data2), np.array(data3), np.array(data4)
        else:
            data = []
            for frame in range(start, stop):
                data.append(self.frames[frame].grab_data(dataname))
            if flatten:
                data = [item for sublist in data for item in sublist]
            if average:
                data = np.mean(data, axis=0)
            return np.array(data)

    def Histogram_Data(self, dataname, molec=2, bins=50, range=(0,100), occupancy=False):
        """ Generates a histogram of the data for the entire simulation

        Args:
            dataname (str): The name of the data to be analyzed.
            molec (int, optional): The molecule number to be analyzed. Defaults to 2.
            bins (int, optional): The number of bins to be used in the histogram. Defaults to 50.
            range (tuple, optional): The range of the histogram. Defaults to (0,100).
            occupancy (bool, optional): Whether or not the data is occupancy data. Defaults to False.

        Returns:
            np.array: The histogram data.
            np.array: The bin edges.

        """
        hist, bin_edges = [], []
        ct = 0
        for frame in self.frames:
            ht, bt = None, None
            if ct == 0:
                ht, bt  = frame.histogram_frame(dataname, molec=molec, bins=bins, range=range, 
                        bin_edges = None, occupancy=occupancy)
                bin_edges = bt
            else:
                ht, bt = frame.histogram_frame(dataname, molec=molec, bins=bins, range=range, 
                        bin_edges = bin_edges, occupancy=occupancy)
            hist.append(ht)
            
            ct = ct + 1
        final_hist = np.sum(np.array(hist), axis=0)
        return final_hist, bin_edges

if __name__ == "__main__":
    data = calculation_data(dirframes=501, dirstart=1, dirstop=100, ncomps=[180,20], compno=[2,6])
    
    pickle.dump(data, open('calc_data.pckl', 'wb'))

    



