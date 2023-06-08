#!/usr/bin/env python
import numpy as np
import pickle

def Read_Hyd_File(data, filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        if (len(lines) > 10):
            for line in lines[3:-1]:
                data.append(line.strip().split())
    return data

def Read_Hbond_File(data, filename, nframes=501):
    mol, h1, h2 = np.genfromtxt(filename, usecols=(0,1,2), unpack=True)
    if data == None:
        data = {}
        data['mol'] = mol.reshape(nframes,-1)
        data['h1'] = h1.reshape(nframes, -1)
        data['h2'] = h2.reshape(nframes, -1)
    else:
        data['mol'] = np.concatenate((data['mol'], mol.reshape(nframes,-1)))
        data['h1'] = np.concatenate((data['h1'], h1.reshape(nframes, -1)))
        data['h2'] = np.concatenate((data['h2'], h2.reshape(nframes, -1)))
    return data



def Pull_Hydration_Data(file_start=1, file_stop=100, fail_okay=False, outname="hydration_data"):
    """This is a python function for pulling the data from a hydration calculation

    Args:
        
    """
    data = []
    for hyd_dir in range(file_start, file_stop):
        try:
            data = Read_Hyd_File(data, "%d/hydration_shell.log"%hyd_dir)
        except:
            if fail_okay:
                print("%d not found, continuing" % hyd_dir)
                continue
            else:
                print("Error: File not found %d" % hyd_dir)
                exit(1)
    final_data = np.asarray(data,dtype=int)
    pickle.dump(final_data, open("%s.pckl"%outname,'wb'))

def Pull_Hbond_Data(file_start=1, file_stop=100, fail_okay=False, outname="hbond_data"):
    data = None
    for hb_dir in range(file_start, file_stop):
        print(hb_dir)
        data = Read_Hbond_File(data, "%d/all_hydrogen_bonds.dat"% hb_dir)
        """
        try: 
            data = Read_Hyd_File(data, "%d/all_hydrogen_bonds.dat"% hb_dir)
        except:
            if fail_okay:
                print("%d not found, continuing" % hb_dir)
            else:
                print("Error: File not found %d" % hb_dir)
                exit(1)
        """
    pickle.dump(data, open("%s.pckl"%outname, 'wb'))

if __name__ == "__main__":
    import argparse
    
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='Choose the options for the analysis setup')

    # Add an argument to the parser
    parser.add_argument('-fstart', default=1, type=int, help='Starting index of calculation')
    parser.add_argument('-fstop', default=100, type=int, help='Stopping index of calculation')
    parser.add_argument('-hydout', default='5_hyd', type=str, help='String to name file with')
    parser.add_argument('-hbout', default='5_hb', type=str, help='String to name file with')
    # Parse the command-line arguments
    args = parser.parse_args()

    Pull_Hydration_Data(file_start=args.fstart, file_stop=args.fstop, fail_okay=True, outname=args.hydout)
    Pull_Hbond_Data(file_start=args.fstart, file_stop=args.fstop, fail_okay=True, outname=args.hbout)
