import MDAnalysis as mda
import pandas as pd
import numpy as np
import math, sys
import argparse
from MDAnalysis.coordinates.memory import MemoryReader
from tqdm import tqdm
from MDAnalysis.analysis import dihedrals
import warnings
from MDAnalysis.analysis import distances

# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

def distance_by_chain(sel1,sel2,chainName):
    distance_46 = distances.distance_array(sel1.positions, # reference
                                    sel2.positions, # configuration
                                    box=u.dimensions)

    dict_distance_46 = {'frame':ts.frame,
                        'chain':chainName,
                        'distance':np.round(distance_46[0],2)
                        }

    df_distance_46 = pd.DataFrame(dict_distance_46)
    return df_distance_46



# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, PDB, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--begin', type=int, default=1,
                    help='Starting Frame. Default = 1 FRAME 1')

parser.add_argument('--end', type=int, default=-1,
                    help='Ending Frame. Default = -1 ALL FRAME')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate for Radial Distribution Function calculations. Default = 1 frames')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Add a system name to output file')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
systemname = args.system

print("LOAD DATA ...")
u = mda.Universe(top_file,traj_file)

L543_A = u.select_atoms("segid PROA and resid 543 and name CG",updating=True)
I637_A = u.select_atoms("segid PROA and resid 637 and name CD",updating=True)
L543_B = u.select_atoms("segid PROB and resid 543 and name CG",updating=True)
I637_B = u.select_atoms("segid PROB and resid 637 and name CD",updating=True)

# print(L543_A)
df = pd.DataFrame()

for ts in tqdm(u.trajectory[::traj_skip]):
    dist_A = distance_by_chain(L543_A,I637_A,"A")
    dist_B = distance_by_chain(L543_B,I637_B,"B")

    df = pd.concat([df,dist_A,dist_B])
    del dist_A,dist_B

# print(df)
df.to_csv("TM46_PROFILE_%s.dat"%(systemname),sep='\t',index=False)
print("\n==> DATA saved to TM46_PROFILE_%s.dat" %(systemname))