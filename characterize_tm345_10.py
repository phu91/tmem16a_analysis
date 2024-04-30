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

def distance_by_chain(frameNow,sel1,sel2,chainName):
    oneDistance = distances.distance_array(sel1.positions, # reference
                                    sel2.positions, # configuration
                                    box=u.dimensions)
    dict_oneDistance = {'frame':frameNow,
                        'resid':sel1.resids,
                        'chain':chainName,
                        'distance':np.round(oneDistance[0],2)
                        }

    df_oneDistance = pd.DataFrame(dict_oneDistance)
    return df_oneDistance



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

R451_A = u.select_atoms("segid PROA and resid 451 and name CZ",updating=True)
K567_A = u.select_atoms("segid PROA and resid 567 and name NZ",updating=True)
K584_A = u.select_atoms("segid PROA and resid 584 and name NZ",updating=True)
K641_A = u.select_atoms("segid PROA and resid 641 and name NZ",updating=True)
K887_A = u.select_atoms("segid PROA and resid 887 and name NZ",updating=True)

R451_B = u.select_atoms("segid PROB and resid 451 and name CZ",updating=True)
K567_B = u.select_atoms("segid PROB and resid 567 and name NZ",updating=True)
K584_B = u.select_atoms("segid PROB and resid 584 and name NZ",updating=True)
K641_B = u.select_atoms("segid PROB and resid 641 and name NZ",updating=True)
K887_B = u.select_atoms("segid PROB and resid 887 and name NZ",updating=True)

df = pd.DataFrame()

for ts in tqdm(u.trajectory[::traj_skip]):
    dist_451_A_887_B = distance_by_chain(ts.frame,R451_A,K887_B,"A")
    dist_567_A_887_B = distance_by_chain(ts.frame,K567_A,K887_B,"A")
    dist_451_B_887_A = distance_by_chain(ts.frame,R451_B,K887_A,"B")
    dist_567_B_887_A = distance_by_chain(ts.frame,K567_B,K887_A,"B")
    dist_564_641_A = distance_by_chain(ts.frame,K584_A,K641_A)
    dist_564_641_B = distance_by_chain(ts.frame,K584_B,K641_B)

    df = pd.concat([df,dist_451_A_887_B,dist_567_A_887_B,dist_451_B_887_A,dist_567_B_887_A,dist_564_641_A,dist_564_641_B])
    del dist_451_A_887_B,dist_567_A_887_B,dist_451_B_887_A,dist_567_B_887_A,dist_564_641_A,dist_564_641_B

# print(df)
with open("TM345_10_PROFILE_%s.dat"%(systemname),'w+') as ofile:
    ofile.write("# DISTANCE BETWEEN R451(A_LOOP23),K567(A_LOOP45) to K887(B_THIRD_CA), K584 to K641\n")

df.to_csv("TM345_10_PROFILE_%s.dat"%(systemname),sep='\t',index=False,mode='a')
print("\n==> DATA saved to TM345_10_PROFILE_%s.dat" %(systemname))