import MDAnalysis as mda
import pandas as pd
import numpy as np
import math, sys
import argparse
from tqdm import tqdm
import warnings

# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

def find_pip2(sel_res_name,sel_res_id):
    PIP2_NEAR_A = u.select_atoms("(around %s (resname %s and resid %s and segid PROA)) and (resname PLPI24 and type PC)"%(cutoff,sel_res_name,sel_res_id),updating=True)
    PIP2_NEAR_B = u.select_atoms("(around %s (resname %s and resid %s and segid PROB)) and (resname PLPI24 and type PC)"%(cutoff,sel_res_name,sel_res_id),updating=True)

    PIP2_NEAR_ID_A = list(set(PIP2_NEAR_A.resids))
    PIP2_NEAR_ID_B = list(set(PIP2_NEAR_B.resids))

    PIP2_NEAR_COUNT_A = len(PIP2_NEAR_ID_A)
    PIP2_NEAR_COUNT_B = len(PIP2_NEAR_ID_B)
    return(PIP2_NEAR_COUNT_A,PIP2_NEAR_COUNT_B)
# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, PDB, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--cutoff', type=float, default=5,
                    help='Cutoff. Default = 5 A')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping frame. Default = 1 frames')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Add a system name to output file')

parser.add_argument('--resname', default='both',nargs='*',
                    help='Select LYS or ARG. Default: BOTH')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
systemname = args.system
cutoff = args.cutoff
select_residue =  args.resname[0].upper()

print("\nLOAD DATA ...\n")
print("Selected RESIDUE: %s\n"%(select_residue))

u = mda.Universe(top_file,traj_file)

if select_residue=='B':
    SELECTED_RESIDUE = u.select_atoms("(resname LYS ARG) and name CA")

else:
    SELECTED_RESIDUE = u.select_atoms("resname %s and name CA"%(select_residue))

SELECTED_RESIDUE_ID = SELECTED_RESIDUE.resids
SELECTED_RESIDUE_NAME = SELECTED_RESIDUE.resnames
# print(SELECTED_RESIDUE_A_ID)
# for i in range(len(SELECTED_RESIDUE_ID_A)):
#     print(SELECTED_RESIDUE_NAME_A[i],SELECTED_RESIDUE_ID_A[i],SELECTED_RESIDUE_NAME_B[i],SELECTED_RESIDUE_ID_B[i])
with open("PIP2_COUNT_CLUSTER_%s_%s.dat"%(systemname,select_residue),'w+') as ofile:
    for ts in tqdm(u.trajectory[::traj_skip],desc='FRAMES'):
        for i in tqdm(range(len(SELECTED_RESIDUE_ID)),desc='Loop through each residue'):
            ofile.write("%s\t%s\t%s\t%s\t%s\n"%(ts.frame,"A",SELECTED_RESIDUE_NAME[i],SELECTED_RESIDUE_ID[i],find_pip2(SELECTED_RESIDUE_NAME[i],SELECTED_RESIDUE_ID[i])[0]))
            ofile.flush()
            ofile.write("%s\t%s\t%s\t%s\t%s\n"%(ts.frame,"B",SELECTED_RESIDUE_NAME[i],SELECTED_RESIDUE_ID[i],find_pip2(SELECTED_RESIDUE_NAME[i],SELECTED_RESIDUE_ID[i])[1]))
            ofile.flush()
