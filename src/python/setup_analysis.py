import numpy as np

def Read_Hyd_File(filename):
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines[2:-1]:
            data.append(line.strip().split())
    return np.as_array(data,dtype=float)



def Pull_Hydration_Data(components=["STYRR","POPC","CLA","SOD","TIP3","LAUR"], to_match=["C","P","N","O"],
                        file_start=1, file_stop=100, fail_okay=False):
    """This is a python function for pulling the data from a hydration calculation

    Args:
        
    """

    for hyd_dir in range(file_start, file_stop):
        Read_Hyd_File("")
