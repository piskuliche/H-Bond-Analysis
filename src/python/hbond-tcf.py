#!/usr/bin/env python

def User_Input():
    import argparse
    parser = argparse.ArgumentParser(description='Code to calculate hbond tcf')
    parser.add_argument('-ncorr',     default=5000,   type=int,   help = 'Correlation length in steps')
    parser.add_argument('-nsep',      default=10,    type=int,   help = 'Separation of origins in steps')
    parser.add_argument('-dt',        default=0.002, type=float, help = 'Simulation timestep in ps')
    parser.add_argument('-dump_freq', default=20,    type=int,   help = 'Simulation dump frequency in steps')
    parser.add_argument('-ndirs',     default=50,    type=int,   help = 'Number of subdirectories')
    parser.add_argument('-fr_per_dir',default=1000,  type=int,   help = 'Number of frames per directory')
    parser.add_argument('-wat_index', default=4096, type=int, help='Starting index for water atoms')
    args = parser.parse_args()

    input_args={}
    input_args['ncorr']       = args.ncorr
    input_args['nsep']        = args.nsep
    input_args['dt']          = args.dt
    input_args['dump_freq']   = args.dump_freq
    input_args['ndirs']       = args.ndirs
    input_args['fpd']         = args.fr_per_dir
    input_args['nframes']     = input_args['ndirs']*input_args['fpd']
    input_args['ntos']        = int((input_args['nframes']-input_args['ncorr'])/input_args['nsep'])
    input_args['wat_index']   = args.wat_index
    return input_args

def Data_Map(input_args):
    data_map=[]
    for i in range(input_args['ndirs']):
        for j in range(input_args['fpd']):
            data_map.append(i+1)
    print(data_map[999:1001])
    return data_map

def Read_Data(input_args,num_dir,hbonds):
    print("Reading directory %d" % num_dir,flush=True)
    hbonds[num_dir] = pd.read_csv("%d/hbond-list.dat"%num_dir, header=None).values
    hbonds[num_dir]=hbonds[num_dir].reshape((input_args['fpd'],-1))
    return hbonds


class All_TCFs:
    def __init__(self,name):
        self.tcfs = {}
    def init_component(self,component):
        self.tcfs[component] = []
        

def TCF_Loop(input_args, hbonds,data_map):
    print("There will be %d origins" % input_args['ntos'])
    data_saved = [0,1]
    tcf = []
    tcf_w2w=[]
    tcf_w2a=[]
    tcf_a2a=[]
    tcf_a2w=[]
    # Outer Loop
    for outer in range(input_args['ntos']):
        t_out = outer*input_args['nsep']
        # Data Management this makes sure that we aren't saving
        # Tons of data all at once.
        while data_map[t_out + input_args['ncorr']] not in data_saved:
            add = max(data_saved)+1
            hbonds=Read_Data(input_args,add,hbonds)
            data_saved.append(add)
        # Clear out data from previous iterations
        hbonds[data_map[t_out]-1] = None

        # Note that the general scheme here is to split into fileblock, and then loc within
        # that block (marked by loc)
        out_loc = t_out%input_args['fpd']
        # Grab hydrogen bonds (shape is # ohs)
        hb_out    = hbonds[data_map[t_out]][out_loc]
        # Initialize a mask that is fully False
        prev_mask = [False]*len(hb_out)
        # Set up origin condition
        out_non_zero = hb_out != 0
        # Calculate number of hydrogen bonds at time t_out
        num_oh_hbond_out = np.sum(out_non_zero*1)
        tcf_temp = []
        # subdivide tcfs into components
        tcf_w2w_temp = []
        tcf_w2a_temp = []
        tcf_a2a_temp = []
        tcf_a2w_temp = []
        # Starting on waters and acceptors
        nz_out_w = hb_out >= input_args['wat_index']
        nz_out_a = hb_out <  input_args['wat_index']
        # Masks for components
        w2w_mask = [False]*len(hb_out)
        w2a_mask = [False]*len(hb_out)
        a2a_mask = [False]*len(hb_out)
        a2w_mask = [False]*len(hb_out)
        # Inner Loop
        for inner in range(input_args['ncorr']):
            t_in = t_out + inner
            loc = t_in%input_args['fpd']
            # Main Calc ---
            # Grab hydrogen bonds (shape is # ohs)
            hb_in = hbonds[data_map[t_in]][loc]
            # Setup masks
            diff         = hb_out != hb_in
            in_non_zero  = hb_in  != 0
            # Make a new mask that includes all criteria, but accounts for things that have 
            # already switched before this (the | (or) criteria)
            eligible    = np.invert(prev_mask)
            prev_mask   = (diff*in_non_zero*out_non_zero)|prev_mask
            # Append to tcf
            tcf_temp.append(np.sum(prev_mask*1))
            # End Main Calc ---

            # Component Calc ---
            is_wat = hb_in >= input_args['wat_index']
            is_acc = (hb_in < input_args['wat_index'])&in_non_zero
            # Update Masks
            w2w_mask = (diff*nz_out_w*is_wat*eligible)|w2w_mask
            w2a_mask = (diff*nz_out_w*is_acc*eligible)|w2a_mask
            a2a_mask = (diff*nz_out_a*is_acc*eligible)|a2a_mask
            a2w_mask = (diff*nz_out_a*is_wat*eligible)|a2w_mask
            # Append to TCF
            tcf_w2w_temp.append(np.sum(w2w_mask*1))
            tcf_w2a_temp.append(np.sum(w2a_mask*1))
            tcf_a2a_temp.append(np.sum(a2a_mask*1))
            tcf_a2w_temp.append(np.sum(a2w_mask*1))
            # End Component Calc ---
        # Append TCFs
        tcf.append(tcf_temp/num_oh_hbond_out)
        tcf_w2w.append(tcf_w2w_temp)
        tcf_w2a.append(tcf_w2a_temp)
        tcf_a2a.append(tcf_a2a_temp)
        tcf_a2w.append(tcf_a2w_temp)
        num_oh_rem = num_oh_hbond_out-tcf_w2w_temp[-1]-tcf_w2a_temp[-1]-tcf_a2a_temp[-1]-tcf_a2w_temp[-1]
        print("Remaining %d" % num_oh_rem,flush=True)
    final_tcf = np.average(tcf,axis=0)
    final_w2w = np.average(tcf_w2w,axis=0)
    final_w2a = np.average(tcf_w2a,axis=0)
    final_a2a = np.average(tcf_a2a,axis=0)
    final_a2w = np.average(tcf_a2w,axis=0)
    times     = np.arange(input_args['ncorr'])*input_args['dump_freq']*input_args['dt']
    np.savetxt('final_tcf.dat',np.c_[times,final_tcf,final_w2w,final_w2a,final_a2a,final_a2w])

    return



if __name__ == "__main__":
    import numpy  as np
    import pandas as pd
    # Read Input
    input_args = User_Input()

    # Set up initial read
    data_map = Data_Map(input_args)
    hbonds={0:None}
    hbonds = Read_Data(input_args,1,hbonds)

    # Calculate TCF
    TCF_Loop(input_args, hbonds, data_map)

    





    

    

