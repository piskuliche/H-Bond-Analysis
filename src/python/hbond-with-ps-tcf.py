#!/usr/bin/env python
import numpy as np

def User_Input():
    import argparse
    parser = argparse.ArgumentParser(description='Code to calculate hbond tcf')
    parser.add_argument('-ncorr',     default=2000,  type=int,   help = 'Correlation length in steps')
    parser.add_argument('-nsep',      default=100,    type=int,   help = 'Separation of origins in steps')
    parser.add_argument('-dt',        default=0.002, type=float, help = 'Simulation timestep in ps')
    parser.add_argument('-dump_freq', default=20,    type=int,   help = 'Simulation dump frequency in steps')
    parser.add_argument('-ndirs',     default=50,    type=int,   help = 'Number of subdirectories')
    parser.add_argument('-fr_per_dir',default=1000,  type=int,   help = 'Number of frames per directory')
    parser.add_argument('-wat_index', default=4096,  type=int,   help = 'Starting index for water atoms')
    parser.add_argument('-num_dir',   default=1,     type=int,   help = 'Directory for read')
    parser.add_argument('-step',      default=1,     type=int,   help = '[0] read [1] run')
    parser.add_argument('-ohyd',      default=1,     type=int,   help = '[0] no ps [1] ps')
    parser.add_argument('-asplit',    default=1,     type=int,   help = '[0] no split [1] split')
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
    input_args['num_dir']     = args.num_dir
    input_args['step']        = args.step
    input_args['ohyd']        = args.ohyd
    input_args['asplit']      = args.asplit
    return input_args

def Data_Map(input_args):
    data_map=[]
    for i in range(input_args['ndirs']):
        for j in range(input_args['fpd']):
            data_map.append(i+1)
    return data_map


class Hbond_TCFS:
    def __init__(self,names):
        self.tcfs = {}
        self.final = {}
        self.mask = {}
        self.temp = {}
        self.eligible = []
        for name in names:
            self.init_component(name)
    def init_component(self,component):
        self.tcfs[component] = []
        self.mask[component] = []
    def new_iter(self,component,n_hb):
        self.temp[component] = []
        self.mask[component] = [False]*n_hb
    def do_comparison(self,component, rules):
        if component == "total":
            self.eligible = np.invert(self.mask["total"])
        rule_prod = self.eligible
        for rule in rules:
            rule_prod = np.multiply(rule_prod,rule)
        self.mask[component] = rule_prod|self.mask[component]
        self.temp[component].append(np.sum(self.mask[component]*1))
    def finish_iter(self,component,num_hb_poss = 1):
        self.tcfs[component].append(np.divide(self.temp[component],num_hb_poss))
    def write_tcf(self,component,nc,df,dt):
        self.times = np.arange(nc)*df*dt
        self.final[component]=np.average(self.tcfs[component],axis=0)
        np.savetxt("%s_hbond_tcf.dat"%component,np.c_[self.times,self.final[component]])


def TCF_Loop(input_args, hbonds, hydrat, data_map):
    print("There will be %d origins" % input_args['ntos'])
    data_saved = [0,1]
    tcf_names = ["total","w2w","w2a","a2a","a2w", "w2w_hyd", "w2w_not","w2a_hyd", "w2a_not"]
    if input_args['asplit'] == 1:
        tcf_names = np.concatenate((tcf_names,['a2a_g', 'a2a_p', 'a2w_g', 'a2w_p']))
    if input_args['ohyd'] == 1:
        tcf_names = np.concatenate((tcf_names,["w2w_ps","w2w_nops","w2a_ps", "w2a_nops"]))
        tcf_names = np.concatenate((tcf_names,["a2a_ps","a2a_nops","a2w_ps", "a2w_nops"]))
        if input_args['asplit'] == 1:
            tcf_names = np.concatenate((tcf_names,["a2a_g_ps","a2a_g_nops", "a2w_g_ps","a2w_g_nops"]))
            tcf_names = np.concatenate((tcf_names,["a2a_p_ps","a2a_p_nops", "a2w_p_ps","a2w_p_nops"]))
    all_tcfs = Hbond_TCFS(tcf_names)

    # <- This next section pulls in the extra data for the extra options, if called
    orig_hyd = []
    hb_g, hb_p = [],[]

    if input_args['asplit'] == 1:
        hb_g = pickle.load(open('hb_g.pckl','rb'))
        hb_p = pickle.load(open('hb_p.pckl','rb'))

    if input_args['ohyd'] == 1:
        orig_hyd = pickle.load(open('ps_and_hyd.pckl','rb'))
        # This checks that the origin hydration is the same shape
        if np.shape(orig_hyd)[0] != int(input_args['nframes']/input_args['nsep']):
            print( np.shape(orig_hyd)[0],int(input_args['nframes']/input_args['nsep']))
            exit("Error: Shapes are not correct!")

    # <- Begin the Outer Loop
    for outer in range(input_args['ntos']):
        t_out = outer*input_args['nsep']
        # Data Management this makes sure that we aren't saving
        # Tons of data all at once.
        while data_map[t_out + input_args['ncorr']] not in data_saved:
            add = max(data_saved)+1
            hbonds=Read_Data_pckl(add,hbonds)
            hydrat=Read_Hydr_pckl(add,hydrat)
            data_saved.append(add)
        # Clear out data from previous iterations
        hbonds[data_map[t_out]-1] = None
        hydrat[data_map[t_out]-1] = None

        # Note that the general scheme here is to split into fileblock, and then loc within
        # that block (marked by loc)
        out_loc = t_out%input_args['fpd']
        # Grab hydrogen bonds (shape is # ohs)
        hb_out    = hbonds[data_map[t_out]][out_loc]
        hyd_out   = hydrat[data_map[t_out]][out_loc]
        # Set up origin condition
        out_non_zero = hb_out != 0
        # Check Hydration
        in_hyd       = hyd_out != 0
        no_hyd       = np.logical_not(in_hyd)
        # Check whether o or g
        ys_g, ys_p = [],[]
        if input_args['asplit'] == 1:
            ys_g = hb_g[outer] != 0
            ys_p = hb_p[outer] != 0
        # Check whether ps is hydrating that lipid
        ys_ps_hyd, no_ps_hyd = [],[]
        if input_args['ohyd'] == 1:
            ys_ps_hyd   = orig_hyd[outer] != 0
            no_ps_hyd   = np.logical_not(ys_ps_hyd)
        # Calculate number of hydrogen bonds at time t_out
        num_oh_hbond_out = np.sum(out_non_zero*1)
        # subdivide tcfs into components
        # Starting on waters and acceptors
        nz_out_w = hb_out >= input_args['wat_index']
        nz_out_a = (hb_out <  input_args['wat_index'])&out_non_zero
        # Masks for components
        for tcf in tcf_names: 
            all_tcfs.new_iter(tcf,len(hb_out))
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
            # End Main calc
            # Component Calc ---
            is_wat = hb_in >= input_args['wat_index']
            is_acc = (hb_in < input_args['wat_index'])&in_non_zero
            # Append to TCF and update masks
            all_tcfs.do_comparison("total",[diff,in_non_zero,out_non_zero])
            all_tcfs.do_comparison("w2w",[diff,nz_out_w,is_wat])
            all_tcfs.do_comparison("w2a",[diff,nz_out_w,is_acc])
            all_tcfs.do_comparison("a2a",[diff,nz_out_a,is_acc])
            all_tcfs.do_comparison("a2w",[diff,nz_out_a,is_wat])
            all_tcfs.do_comparison("w2w_hyd",[diff,nz_out_w,is_wat,in_hyd])
            all_tcfs.do_comparison("w2a_hyd",[diff,nz_out_w,is_acc,in_hyd])
            all_tcfs.do_comparison("w2w_not",[diff,nz_out_w,is_wat,no_hyd])
            all_tcfs.do_comparison("w2a_not",[diff,nz_out_w,is_acc,no_hyd])
            if input_args["ohyd"] == 1:
                all_tcfs.do_comparison("w2w_ps",[diff,nz_out_w,is_wat,ys_ps_hyd,in_hyd])
                all_tcfs.do_comparison("w2w_nops",[diff,nz_out_w,is_wat,no_ps_hyd,in_hyd])
                all_tcfs.do_comparison("w2a_ps",[diff,nz_out_w,is_acc,ys_ps_hyd,in_hyd])
                all_tcfs.do_comparison("w2a_nops",[diff,nz_out_w,is_acc,no_ps_hyd,in_hyd])
                all_tcfs.do_comparison("a2a_ps",[diff,nz_out_a,is_acc,ys_ps_hyd])
                all_tcfs.do_comparison("a2a_nops",[diff,nz_out_a,is_acc,no_ps_hyd])
                all_tcfs.do_comparison("a2w_ps",[diff,nz_out_a,is_wat,ys_ps_hyd])
                all_tcfs.do_comparison("a2w_nops",[diff,nz_out_a,is_wat,no_ps_hyd])
            if input_args['asplit'] == 1:
                all_tcfs.do_comparison("a2a_g",[diff,nz_out_a,is_acc,ys_g])
                all_tcfs.do_comparison("a2a_p",[diff,nz_out_a,is_acc,ys_p])
                all_tcfs.do_comparison("a2w_g",[diff,nz_out_a,is_wat,ys_g])
                all_tcfs.do_comparison("a2w_p",[diff,nz_out_a,is_wat,ys_p])
                if input_args['ohyd'] == 1:
                    all_tcfs.do_comparison("a2a_g_ps", [diff,nz_out_a,is_acc,ys_g,ys_ps_hyd])
                    all_tcfs.do_comparison("a2a_p_ps", [diff,nz_out_a,is_acc,ys_p,ys_ps_hyd])
                    all_tcfs.do_comparison("a2w_g_ps", [diff,nz_out_a,is_wat,ys_g,ys_ps_hyd])
                    all_tcfs.do_comparison("a2w_p_ps", [diff,nz_out_a,is_wat,ys_p,ys_ps_hyd])
                    all_tcfs.do_comparison("a2a_g_nops",[diff,nz_out_a,is_acc,ys_g,no_ps_hyd])
                    all_tcfs.do_comparison("a2a_p_nops",[diff,nz_out_a,is_acc,ys_p,no_ps_hyd])
                    all_tcfs.do_comparison("a2w_g_nops",[diff,nz_out_a,is_wat,ys_g,no_ps_hyd])
                    all_tcfs.do_comparison("a2w_p_nops",[diff,nz_out_a,is_wat,ys_p,no_ps_hyd])
            # End Component Calc ---
        # Append TCFs
        all_tcfs.finish_iter("total",num_oh_hbond_out)
        for tcf in tcf_names[1:]:
            all_tcfs.finish_iter(tcf)
    pickle.dump(all_tcfs,open('all_tcfs.pckl','wb'))
    for tcf in tcf_names:
        all_tcfs.write_tcf(tcf,input_args["ncorr"],input_args["dump_freq"],input_args["dt"])
    return

def Read_and_Convert(fpd, num_dir):
    hbonds = pd.read_csv("%d/hbond-list.dat"%num_dir, header=None).values
    hbonds = hbonds.reshape((fpd,-1))
    pickle.dump(hbonds,open("%d/hbonds.pckl"%num_dir,'wb'))
    return

def Read_and_Convert_Hydration(fpd,num_dir):
    hydrate = pd.read_csv("%d/hydration-list.dat"%num_dir, header=None,delim_whitespace=True).values.T[0]
    hydrate = np.array(list(map(int,hydrate)))
    hydrate = hydrate.reshape((fpd,-1))
    pickle.dump(hydrate,open("%d/hydration.pckl"%num_dir,'wb'))
    return

def Read_Data_pckl(num_dir,hbonds):
    print("Reading directory %d" % num_dir,flush=True)
    hbonds[num_dir] = pickle.load(open('%d/hbonds.pckl'%num_dir,'rb'))
    return hbonds

def Read_Hydr_pckl(num_dir,hydrat):
    print("Reading directory %d" % num_dir,flush=True)
    hydrat[num_dir] = pickle.load(open('%d/hydration.pckl'%num_dir,'rb'))
    return hydrat


if __name__ == "__main__":
    import numpy  as np
    import pandas as pd
    import pickle
    # Read Input
    input_args = User_Input()

    if input_args["step"] == 0:
        Read_and_Convert(input_args["fpd"],input_args["num_dir"])
        Read_and_Convert_Hydration(input_args["fpd"],input_args["num_dir"])
    else:
        # Set up initial read
        data_map = Data_Map(input_args)
        hbonds= {0:None}
        hydrat= {0:None}
        hbonds = Read_Data_pckl(1,hbonds)
        hydrat = Read_Hydr_pckl(1,hydrat)

        # Calculate TCF
        TCF_Loop(input_args, hbonds, hydrat, data_map)

    





    

    

