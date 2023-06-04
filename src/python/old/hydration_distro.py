#!/usr/bin/env python3
import numpy as np

def grab_hydration(ndirs, nrows):
    """Grabs the hydration data from each of n directories

    Grabs data from n/hydration_shell.log where n is an integer from 1 to n, where there
    are n total files. Each file has 9 columns, the first being the timestep, 
    and the rest are the hydration data.

    Args:
        ndirs (int): Number of directories to grab data from
        nrows (int): Number of rows in each file

    Returns:
        np.array: Array of shape (8, nrows*ndirs) containing the data
    """
    data = np.zeros((8, nrows*ndirs))
    for i in range(ndirs):
        vals = np.genfromtxt(f'{i+1}/hydration_shell.log', skip_header=1, skip_footer=1, max_rows=nrows)
        data[:, i*nrows:(i+1)*nrows] = vals.T
    return data



if __name__ == "__main__":
    import argparse
    import pickle

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', default=20, type=int, help='Number of directories')
    Iargs = parser.parse_args()

    data = grab_hydration(Iargs.n, 1000)

    pickle.dump(data, open('hydration_counts.pckl', 'wb'))
