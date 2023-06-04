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
    isHeavy=[]
    LaurisHeavy=[]
    for atom in u.atoms:
        if ("O" in atom.name) and atom.resname == "DOPC":
            isAcc.append(1)
        elif atom.resname == "DOPC":
            isAcc.append(0)
    for atom in u.atoms:
        if ("C" in atom.name) and atom.resname == "DOPC":
            isHeavy.append(1)
        elif ("O" in atom.name) and atom.resname == "DOPC":
            isHeavy.append(2)
        elif ("P" in atom.name) and atom.resname == "DOPC":
            isHeavy.append(3)
        elif ("N" in atom.name) and atom.resname == "DOPC":
            isHeavy.append(4)
        elif atom.resname == "DOPC":
            isHeavy.append(0)
    for atom in u.atoms:
        if ("C" in atom.name) and atom.resname == "LAUR":
            LaurisHeavy.append(5)
        elif ("O" in atom.name) and atom.resname == "LAUR":
            LaurisHeavy.append(6)
        elif ("P" in atom.name) and atom.resname == "LAUR":
            LaurisHeavy.append(7)
        elif ("N" in atom.name) and atom.resname == "LAUR":
            LaurisHeavy.append(8)
        elif atom.resname == "LAUR":
            LaurisHeavy.append(0)

    numacc = np.sum(isAcc)
    f=open('is_acc.txt','w')
    f.write("%d %d\n" % (numacc,len(isAcc)))
    for i in range(len(isAcc)):
        f.write("%d\n" % isAcc[i])
    f.close()
    f=open('is_heavy.txt','w')
    f.write("%d\n" % len(isHeavy))
    for i in range(len(isHeavy)):
        f.write("%d\n" % isHeavy[i])
    f.close()
    f=open('laur_is_heavy.txt','w')
    f.write("%d\n" % len(LaurisHeavy))
    for i in range(len(LaurisHeavy)):
        f.write("%d\n" % LaurisHeavy[i])
    f.close()
    f=open('heavy.counts','w')
    isHeavy=np.array(isHeavy)
    Ccount = np.sum((isHeavy==1)*1)
    Ocount = np.sum((isHeavy==2)*1)
    Pcount = np.sum((isHeavy==3)*1)
    Ncount = np.sum((isHeavy==4)*1)
    f.write("%d %d %d %d\n" % (Ccount,Ocount,Pcount,Ncount))
    f.close()
    f=open('laur.counts','w')
    LaurisHeavy=np.array(LaurisHeavy)
    Ccount = np.sum((LaurisHeavy==5)*1)
    Ocount = np.sum((LaurisHeavy==6)*1)
    Pcount = np.sum((LaurisHeavy==7)*1)
    Ncount = np.sum((LaurisHeavy==8)*1)
    f.write("%d %d %d %d\n" % (Ccount,Ocount,Pcount,Ncount))
    f.close()


