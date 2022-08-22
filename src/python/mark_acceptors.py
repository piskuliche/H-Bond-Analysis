#!/usr/bin/env python
"""
This is a short python code that is part of the repository git@github.com:piskuliche/Contact-Analysis.git

All it does is takes a gromacs-style TPR file (as an argument) and then writes a file called "is_atom.txt" in the working directory to mark whether that atom is in group 1, group 2, or neither. Group 1 is the analysis group.

Note - the if statement MUST be modified to set criteria.

Questions? Contact piskulichz@gmail.com

Copyright Aug 2022, Zeke Piskulich

"""
import MDAnalysis as mda
import sys
import numpy as np
import argparse




if __name__ == "__main__":
    parser=argparse.ArgumentParser(description='Code to mark atoms for contact analysis')
    parser.add_argument("-tpr",default="topol.tpr",type=str,help='Gromacs topology file')
    args = parser.parse_args()
    tprfile = args.tpr
    print(tprfile)

    u = mda.Universe(tprfile)

    isAcc= []
    for atom in u.atoms:
        if ("O" in atom.name) and atom.resname == "DOPC":
            isAcc.append(1)
        elif atom.resname == "DOPC":
            isAcc.append(0)

    numacc = np.sum(isAcc)
    f=open('is_acc.txt','w')
    f.write("%d %d\n" % (numacc,len(isAcc)))
    for i in range(len(isAcc)):
        f.write("%d\n" % isAcc[i])
    f.close()

