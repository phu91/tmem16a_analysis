import MDAnalysis as mda
import pandas as pd
import numpy as np
import math, sys
import argparse
from tqdm import tqdm
import warnings

import matplotlib.pyplot as plt
import seaborn as sns 
from matplotlib import rc

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

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping frame. Default = 1 frames')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Add a system name to output file')

parser.add_argument('--nbins', type=int, default=5,
                    help='Number of bins')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
systemname = args.system
nbins = args.nbins


print("\nLOAD DATA ...\n")

u = mda.Universe(top_file,traj_file)

P_head = u.select_atoms("type PL",updating=True)
n_lipid = 660
n_frame = len(u.trajectory[::traj_skip])
# print(n_frame)

P_DIMENSION_X_LIST = []
P_DIMENSION_Y_LIST = []

for ts in u.trajectory[::traj_skip*50]:  ## TESTING CASE. IN PRODUCTION RUN, traj_skip needs to be timed a higher factor
    P_head = u.select_atoms("type PL",updating=True)
    P_DIMENSION = P_head.positions
    for i in range(n_lipid):
        P_DIMENSION_X_LIST=np.append(P_DIMENSION[i][0],P_DIMENSION_X_LIST)
        P_DIMENSION_Y_LIST=np.append(P_DIMENSION[i][1],P_DIMENSION_Y_LIST)

# print(P_DIMENSION_X_LIST.min(),P_DIMENSION_X_LIST.max())
# print(P_DIMENSION_Y_LIST.min(),P_DIMENSION_Y_LIST.max())

nbins=nbins
GRID_X = np.histogram_bin_edges(P_DIMENSION_X_LIST,bins=nbins)
GRID_Y = np.histogram_bin_edges(P_DIMENSION_Y_LIST,bins=nbins)
# print(GRID_X)
# print(GRID_Y)
THICKNESS_LIST = []

with open("MEMBRANE_THICKNESS_COM_%s.dat"%(systemname),'w+') as comfile:
    comfile.write("#FRAME COM_A_X COM_A_Y COM_A_Z COM_B_X COM_B_Y COM_B_Z\n")
    for ts in tqdm(u.trajectory[::traj_skip],desc='FRAMES'):
        for i in range(nbins):
            for j in range(nbins):
                P_head_selected_upper = u.select_atoms("type PL and resid 1 to 400 and (prop x >= %s and prop x < %s and prop y >= %s and prop y < %s)"%(GRID_X[i],GRID_X[i+1],GRID_Y[j],GRID_Y[j+1]),updating=True)        
                P_head_selected_lower = u.select_atoms("type PL and resid 401 to 660 and (prop x >= %s and prop x < %s and prop y >= %s and prop y < %s)"%(GRID_X[i],GRID_X[i+1],GRID_Y[j],GRID_Y[j+1]),updating=True)        
                P_head_selected_upper_ID = P_head_selected_upper.resids
                P_head_selected_lower_ID = P_head_selected_lower.resids
                # print(P_head_selected_upper_ID,P_head_selected_lower_ID)
                if len(P_head_selected_upper_ID)!=0 and len(P_head_selected_lower_ID)!=0:
                    P_head_selected_positions_upper = P_head_selected_upper.center_of_mass()
                    P_head_selected_positions_lower = P_head_selected_lower.center_of_mass()
                    thickness = P_head_selected_positions_upper[2]-P_head_selected_positions_lower[2]
                    thickness = round(thickness,3)
                else:
                    thickness=0.000   # AREA WHICH HAS NO LIPID
                THICKNESS_LIST=np.append((int(ts.frame),round((GRID_X[i]+GRID_X[i+1])/2,3),round((GRID_Y[j]+GRID_Y[j+1])/2,3),thickness),THICKNESS_LIST)
    # print(THICKNESS_LIST)
        chainA = u.select_atoms("segid PROA and resid 400:700")
        chainB = u.select_atoms("segid PROB and resid 400:700")
        chainA_COM = chainA.center_of_mass()
        chainB_COM = chainB.center_of_mass()
        comfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(ts.frame,*chainA_COM,*chainB_COM))
        comfile.flush()
THICKNESS_LIST = np.reshape(THICKNESS_LIST,(n_frame,nbins**2,4))
# print(THICKNESS_LIST)

THICKNESS_LIST_DF = pd.DataFrame(columns=['FRAME','X','Y','Z_THICK'])
for i in range(n_frame):
    for j in range(nbins**2):
        _ = tuple(THICKNESS_LIST[i][j])
        _DF = pd.DataFrame([_],columns=['FRAME','X','Y','Z_THICK'])
        THICKNESS_LIST_DF = pd.concat([THICKNESS_LIST_DF,_DF])
THICKNESS_LIST_DF = THICKNESS_LIST_DF.sort_values(['FRAME']).reset_index()
THICKNESS_LIST_DF = THICKNESS_LIST_DF[['FRAME','X','Y','Z_THICK']]
# print(THICKNESS_LIST_DF)
with open("MEMBRANE_THICKNESS_%s.dat"%(systemname),'w+') as ofile:
    ofile.write("#BINS: %s\n"%(nbins))
    ofile.write("#FRAME X Y Z_THICKNESS\n")
    ofile.write("\n")

THICKNESS_LIST_DF.to_csv("MEMBRANE_THICKNESS_%s.dat"%(systemname),
                        mode='a',
                        index=False,sep='\t',header=False,
                        chunksize=10000)