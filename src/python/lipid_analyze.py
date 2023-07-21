#!/usr/bin/env python

import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations
import matplotlib.pyplot as plt
import math


"""
**** 
Basic Math Functions
****
"""
def P2(x):
    return (3*np.cos(x)**2. - 1)/2.0

def bond_unit_vec(r1, r2, L):
    dr = r2 - r1
    dr = dr - L*np.rint(dr/L)
    return dr/np.sum(dr**2.)

"""
**** 
Data Storage Class
****
"""

class MolData:
    def __init__(self, molecs, data, reference):
        self.molids = reference
        self.data = np.zeros_like(reference)
        self.add_by_reference(molecs, data, reference)
        return
    def add_by_reference(self, molecs, data, reference):
        mask = np.isin(reference, molecs)
        self.data[mask] = data[np.where(molecs == reference[mask])]
        return
        

            
    
"""
**** 
Analysis Functions
****
"""

def Analyze_by_Residue(u, ts, AtomGroup, function, optional=None):
    residues = np.unique(AtomGroup.resids)
    data = np.zeros(len(residues))
    molids = np.zeros(len(residues))
    for i, residue in enumerate(residues):
        ResGroups = AtomGroup.split("residue")
        data[i] = function(u, ts, AtomGroup, ResGroups[i], optional=optional)
        molids[i] = ResGroups[i].resids[0]
    return molids, data

def Calc_Z_byres(u, ts, AtomGroup, ResAtoms, optional = None):
    if optional != None:
        raise ValueError("Optional must be None for Calc_Z_byres")
    zcom = AtomGroup.center_of_mass()[2]
    zdist = ResAtoms.center_of_mass(wrap=True)[2] - zcom
    return zdist

def Calc_Rg_byres(u, ts, AtomGroup, ResAtoms, optional = None):
    if optional != None:
        raise ValueError("Optional must be None for Calc_Rg_byres")
    return ResAtoms.radius_of_gyration()

def Calc_P2_byres(u, ts, AtomGroup, ResAtoms, optional = [0,0,1]):
    if len(ResAtoms) != 2:
        print("ResAtoms must be only two atoms")

    L = ts.dimensions[:3]
    pos1 = ResAtoms.positions[0]
    pos2 = ResAtoms.positions[1]
    bond_vec = bond_unit_vec(pos1, pos2, L)
    cosangle = np.dot(bond_vec, optional)
    p2_value = P2(cosangle)
    return p2_value

def Calc_CosAngle_byres(u, ts, AtomGroup, ResAtoms, optional = [0,0,1]):
    if len(ResAtoms) != 2:
        print("ResAtoms must be only two atoms")

    L = ts.dimensions[:3]
    pos1 = ResAtoms.positions[0]
    pos2 = ResAtoms.positions[1]
    bond_vec = bond_unit_vec(pos1, pos2, L)
    cosangle = np.dot(bond_vec, optional)
    return cosangle

def Calc_Angle_byres(u, ts, AtomGroup, ResAtoms, optional = [0,0,1]):
    if len(ResAtoms) != 2:
        print("ResAtoms must be only two atoms")

    L = ts.dimensions[:3]
    pos1 = ResAtoms.positions[0]
    pos2 = ResAtoms.positions[1]
    bond_vec = bond_unit_vec(pos1, pos2, L)
    cosangle = np.dot(bond_vec, optional)
    angle = np.arccos(cosangle)
    return angle

def Main_Analysis(output_data, tprfile, trajfile):
    u = mda.Universe(tprfile, trajfile)
    workflow = [transformations.unwrap(u.atoms)]
    u.trajectory.add_transformations(*workflow)
    membrane = u.select_atoms("resname POPC or resname LAUR")
    laur_CN = u.select_atoms("resname LAUR and (type CG2O5 or type NG301)")
    
    for ts in u.trajectory:
        print(ts.time)
        L = ts.dimensions
        zmol, zdat = Analyze_by_Residue(u, ts, membrane, Calc_Z_byres)
        output_data["zdists"].append(MolData(zmol, zdat, zmol))
        rmol, rdat = Analyze_by_Residue(u, ts, membrane, Calc_Rg_byres)
        output_data["rg"].append(MolData(rmol, rdat, zmol))
        pmol, pdat = Analyze_by_Residue(u, ts, laur_CN, Calc_P2_byres, optional=[0,0,1])
        print(np.shape(zdat))
        output_data["tiltp2"].append(MolData(pmol,pdat, zmol))
        cmol, cdat = Analyze_by_Residue(u, ts, laur_CN, Calc_CosAngle_byres, optional=[0,0,1])
        output_data["tiltcosangle"].append(MolData(cmol,cdat, zmol))
        amol, adat = Analyze_by_Residue(u, ts, laur_CN, Calc_Angle_byres, optional=[0,0,1])
        output_data["tiltangle"].append(MolData(amol,adat, zmol))

    return output_data

def Do_Files(start=1,stop=100):
    output_data = {"zdists":[], "rg":[], "tiltangle":[], "tiltcosangle":[], "tiltp2":[]}
    for i in range(start,stop):
        tprfile = "../../tpr/step7_%d.tpr"%i
        xtcfile = "../../xtc/step7_%d.xtc"%i
        output_data = Main_Analysis(output_data, tprfile, xtcfile)

def Do_File(filenumber=1):
    tprfile = "../../tpr/step7_%d.tpr"%filenumber
    xtcfile = "../../xtc/step7_%d.xtc"%filenumber
    output_data = {"zdists":[], "rg":[], "tiltangle":[], "tiltcosangle":[], "tiltp2":[]}
    output_data = Main_Analysis(output_data, tprfile, xtcfile)
    Return_Data(output_data, filenumber)

def Return_Data(output_data, number):
    import pickle
    cnt = 0
    for key in output_data:
        outname = "output_%s_%s.pckl"%(key, number)
        dim2_data, molids = [], []
        for i, frame in enumerate(output_data[key]):
            if i == 0:
                molids = frame.molids
            dim2_data.append(frame.data)
        pickle.dump(dim2_data, open(outname, 'wb'))
        if cnt == 0: pickle.dump(molids, open("molids_%s_%s.pckl"%(key, number), 'wb'))
        cnt = 1
    return

if __name__ == "__main__":
    import argparse
    import sys

    filenumber = int(sys.argv[1])

    Do_File(filenumber)