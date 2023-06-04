#!/usr/bin/env python3
import numpy as np

def Pull_Atom_Names(gro_file, lipid_name):
    """Pulls the atom names from a gro file and returns them as an ordered list.

    Args: 
        gro_file (str): The path to the gro file.
        lipid_name (str): The name of the lipid in the gro file.

    """
    
    atom_names = []
    with open(gro_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if lipid_name in line:
                temporary = line.strip().split()
                if "C" in temporary[1]:
                    atom_names.append(1)
                elif "O" in temporary[1]:
                    atom_names.append(2)
                elif "P" in temporary[1]:
                    atom_names.append(3)
                elif "N" in temporary[1]:
                    atom_names.append(4)
                else:
                    atom_names.append(0)

    with open('is_heavy.txt', 'w') as f:
        f.write('%d\n' % len(atom_names))
        for atom in atom_names:
            f.write('%d\n' % atom)

    return atom_names

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", default='step8_1.gro', type=str, help="The path to the gro file.")
    parser.add_argument("-lipid", default='DOPC', type=str, help="The name of the lipid in the gro file.")
    args = parser.parse_args()

    Pull_Atom_Names(args.g, args.lipid)