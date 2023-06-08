#!/usr/bin/env python
import numpy as np
import MDAnalysis as mda

def Read_GRO_and_Match(grofile, contribute_matches=["POPC","LAUR","TIP3"], to_match=["C","P","N","O"],
                        outfile="map_of_heavy_atoms.info"):
    """This is a python function that reads the gro file and finds matches
    """
    u = mda.Universe(grofile)
    all_atoms = u.select_atoms("all")
    to_match=["P", "O", "N", "C"]
    nu_match=[1,2,3,4]
    atom_labels = []
    per_type_match = {}
    per_type_nomatch = {}
    for key in contribute_matches:
        per_type_match[key] = 0
    for i, atype in enumerate(all_atoms.types):
        found = False
        resname = all_atoms[i].resname
        if all_atoms[i].resname in contribute_matches:
            for j, key in enumerate(to_match):
                if key in atype:
                    per_type_match[resname] += 1
                    atom_labels.append(nu_match[j])
                    found = True
        if found == False:
            atom_labels.append(0)
            if resname not in per_type_nomatch:
                per_type_nomatch[resname] = 0
            per_type_nomatch[resname] += 1

    print("There are %d atoms that match" % np.sum(atom_labels))
    print("*** Residues and How Many Matches ***")
    print(per_type_match)
    print("*** Residues and How Many Non Matches ***")
    print(per_type_nomatch)
    np.savetxt(outfile, np.c_[atom_labels], fmt='%i')

def Read_GRO_and_Match_ACC(grofile, contribute_matches=["POPC","LAUR","TIP3"], to_match=["N","O"], nu_match=[1,1],
                        outfile="map_of_acc_atoms.info"):
    """This is a python function that reads the gro file and finds matches
    """
    u = mda.Universe(grofile)
    all_atoms = u.select_atoms("all")

    atom_labels = []
    per_type_match = {}
    per_type_nomatch = {}
    for key in contribute_matches:
        per_type_match[key] = 0
    for i, atype in enumerate(all_atoms.types):
        found = False
        resname = all_atoms[i].resname
        if all_atoms[i].resname in contribute_matches:
            for j, key in enumerate(to_match):
                if key in atype:
                    per_type_match[resname] += 1
                    atom_labels.append(nu_match[j])
                    found = True
        elif all_atoms[i].resname == "TIP3":
            if "OH2" in atype:
                per_type_match[resname] += 1
                atom_labels.append(2)
                found=True
        if found == False:
            atom_labels.append(0)
            if resname not in per_type_nomatch:
                per_type_nomatch[resname] = 0
            per_type_nomatch[resname] += 1

    print("There are %d atoms that match" % np.sum(atom_labels))
    print("*** Residues and How Many Matches ***")
    print(per_type_match)
    print("*** Residues and How Many Non Matches ***")
    print(per_type_nomatch)
    np.savetxt(outfile, np.c_[atom_labels], fmt='%i')
        



if __name__ == "__main__":
    print("test")
    Read_GRO_and_Match("gro/step7_1.gro")
    Read_GRO_and_Match_ACC("gro/step7_1.gro")
