#!/usr/bin/env python

import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations
import matplotlib.pyplot as plt


def Calculate_Z_Distance(membrane):
    zcom = membrane.center_of_mass(wrap=True)[2]
    zdists = np.abs(membrane.center_of_mass(wrap=True, compound='residues')[:,2] - zcom)
    return zdists

def Calculate_Rg(u, AtomGroup):
    residues = np.unique(AtomGroup.resids)
    rg_vals = []
    for residue in residues:
        molgroup = u.select_atoms("resid %d" % residue)
        rg_vals.append(molgroup.radius_of_gyration())
    print(np.average(rg_vals), np.min(rg_vals), np.max(rg_vals))
    return rg_vals
    

def Main_Analysis(grofile, trajfile):

    u = mda.Universe(grofile, trajfile)
    workflow = [transformations.unwrap(u.atoms)]
    u.trajectory.add_transformations(*workflow)
    membrane = u.select_atoms("resname POPC or resname LAUR")
    popc_molecs = u.select_atoms("resname POPC")
    laur_molecs = u.select_atoms("resname LAUR")
    output_data = {"zdists":[], "rg":[]}
    for ts in u.trajectory:
        output_data["zdists"].append(Calculate_Z_Distance(membrane))
        output_data["rg"].append(Calculate_Rg(u, popc_molecs))
        
    return


if __name__ == "__main__":
    import argparse

    Main_Analysis("../tpr/step7_50.tpr", "../xtc/step7_50.xtc")