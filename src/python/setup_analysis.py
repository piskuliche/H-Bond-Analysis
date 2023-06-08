import numpy as np
import pickle

def Read_Hyd_File(data, filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        if (len(lines) > 10):
            for line in lines[3:-1]:
                data.append(line.strip().split())
    return data



def Pull_Hydration_Data(components=["STYRR","POPC","CLA","SOD","TIP3","LAUR"], to_match=["C","P","N","O"],
                        file_start=1, file_stop=100, fail_okay=False, outname="hydration_data"):
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

if __name__ == "__main__":
    Pull_Hydration_Data(file_start=1, file_stop=80, fail_okay=True, outname="10_hyd")
