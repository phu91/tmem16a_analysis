# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.visualization import streamlines
import matplotlib.pyplot as plt
import math, sys
from MDAnalysis.lib.distances import distance_array
import seaborn as sns 
import argparse
from tqdm import tqdm

def find_chlorides():
    # print("\n===> CUT-OFF distance of CHLORIDES and the protein is %s Angstroms\n"%(cutoff))
    cutoff=3
    traj_skip_fast = 50
    chloride_list_found = []
    # with mda.Writer('CHLORIDES_at_pore_chain_A_%s.dcd'%(systemname), all.n_atoms) as w:
    for ts in tqdm(u.trajectory[::traj_skip_fast]):
        # print(ts.frame)
        chloride_a = u.select_atoms("(around %s (resid 584 641 599)) and resname CLA"%(cutoff),updating=True)
        if len(chloride_a)!=0:
            # print(chloride_a)
            for ion in chloride_a.atoms.indices:
                chloride_list_found.append((ion))
                # w.write(all)

    # with mda.Writer('CHLORIDES_at_pore_chain_B_%s.dcd'%(systemname), all.n_atoms) as w:
        # for ts in tqdm(u.trajectory[::traj_skip_fast],desc='Finding Cl- at chain B'):
        #     # print(ts.frame)
        #     chloride_b = u.select_atoms("(around %s (resid 584 641 599 and segid PROB)) and resname CLA"%(cutoff),updating=True)
        #     if len(chloride_b)!=0:
        #         for ion in chloride_b.atoms.indices:
        #             chloride_list_found.append((ion,"B"))
        #             w.write(all)
    # print("\n===> Trajectories of FOUND Cl- are saved in 'Chlorides_at_pore_chain_A_%s.dcd' and 'Chlorides_at_pore_chain_B_%s.dcd'\n"%(systemname,systemname))
    return chloride_list_found

# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate. Default = 1 frames')

parser.add_argument('--resid', type=int, nargs='+',
                    help='List of selected Chloride')

args = parser.parse_args()

top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
selected_resid=args.resid

u = mda.Universe(top_file,traj_file)
print(u)
step_lys_1 = u.select_atoms("resname LYS and resid 584",updating=True).center_of_mass()[2]
step_lys_2 = u.select_atoms("resname LYS and resid 641",updating=True).center_of_mass()[2]
step_lys_3 = u.select_atoms("resname LYS and resid 599",updating=True).center_of_mass()[2]

if selected_resid is None:
    print("\n===> Activate FIND_CHLORIDE automatically!\n")
    found_chlorides=find_chlorides()
    if len(found_chlorides)==0:
        print("\n===> No Chloride found at the Pore! Please input a specific index for Chloride.\n")
        breakpoint
    else:
        selected_resid=found_chlorides
        with open("HYDRATION_Z_Cl_at_PORE.dat",'w+') as of:
            of.write("#Frame | Cl_ID | Hydration | Zcoordinate | STEP1 | STEP2 |STEP3\n")
            for ts in u.trajectory[::traj_skip]:
                for cl in selected_resid:
                    selected_chloride = u.select_atoms("name CLA and index %s"%(cl),updating=True)
                    hydration_water = u.select_atoms("(around 3.2 (name CLA and index %s)) and name OH2"%(cl),updating=True)
                    print("frame %s ID: %s"%(ts.frame,cl))
                    of.write("%i\t%i\t%i\t%2.2f\t%2.2f\t%2.2f\t%2.2f\n"%(ts.frame,cl,len(hydration_water),*selected_chloride.positions[:,2],step_lys_1,step_lys_2,step_lys_3))
else:
    with open("HYDRATION_Z_Cl_SELECTED.dat",'w+') as of:
        of.write("#Frame | Cl_ID | Hydration | Zcoordinate | STEP1 | STEP2 |STEP3\n")
        for ts in u.trajectory[::traj_skip]:
            for cl in selected_resid:
                selected_chloride = u.select_atoms("name CLA and index %s"%(cl),updating=True)
                hydration_water = u.select_atoms("(around 3.2 (name CLA and index %s)) and name OH2"%(cl),updating=True)
                print("frame %s ID: %s"%(ts.frame,cl))
                print(ts.frame,cl,len(hydration_water),*selected_chloride.positions[:,2],step_lys_1,step_lys_2,step_lys_3)
                # of.write("%i\t%i\t%i\t%2.2f\t%2.2f\t%2.2f\t%2.2f\n"%(ts.frame,cl,len(hydration_water),*selected_chloride.positions[:,2],step_lys_1,step_lys_2,step_lys_3))
