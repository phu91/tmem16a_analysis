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

def extract_frames(startFrame,endFrame):
    protein = u.select_atoms("protein and name CA")
    with mda.Writer("SEGMENT_FRAME_%s_TO_%s.pdb"%(startFrame,endFrame)) as pdb:
        pdb.write(protein)

    with mda.Writer("SEGMENT_FRAME_%s_TO_%s.xtc"%(startFrame,endFrame), protein.n_atoms) as W:
        for ts in u.trajectory[startFrame:endFrame:traj_skip]:
            W.write(protein)

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
traj_begin = args.begin
traj_end = args.end
systemname = args.system

print("LOAD DATA ...")
u = mda.Universe(top_file,traj_file)
n_atom_origin = len(u.atoms)


if traj_end != -1:
    extract_frames(traj_begin,traj_end)
    u = mda.Universe("SEGMENT_FRAME_%s_TO_%s.pdb"%(traj_begin,traj_end),
                     "SEGMENT_FRAME_%s_TO_%s.xtc"%(traj_begin,traj_end),
                     in_memory=True)
    print("\n########################################################")
    print("TOP : SEGMENT_FRAME_%s_TO_%s.pdb"%(traj_begin,traj_end))
    print("TRAJ: SEGMENT_FRAME_%s_TO_%s.xtc"%(traj_begin,traj_end))
    print("NUMBER OF ATOMS (ORIGINAL): %s"%(n_atom_origin))
    print("NUMBER OF ATOMS (CURRENT) : %s"%(len(u.atoms)))
    print("NUMBER OF FRAMES: %s"%(len(u.trajectory)))
    print("########################################################\n")
else:
    print("\n########################################################")
    print("TOP : %s"%(top_file))
    print("TRAJ: %s"%(traj_file))
    print("NUMBER OF ATOMS (ORIGINAL): %s"%(n_atom_origin))
    print("NUMBER OF ATOMS (CURRENT) : %s"%(len(u.atoms)))
    print("NUMBER OF FRAMES: %s"%(len(u.trajectory)))
    print("########################################################\n")
    pass

ref_850_A = u.select_atoms("protein and resid 850 and name CA and segid PROA",updating=True)
ref_850_B = u.select_atoms("protein and resid 850 and name CA and segid PROB",updating=True)
ref_887_A = u.select_atoms("protein and resid 887 and name CA and segid PROA",updating=True)
ref_887_B = u.select_atoms("protein and resid 887 and name CA and segid PROB",updating=True)
res_889_A = u.select_atoms("protein and resid 889 and name CA and segid PROA",updating=True)
res_889_B = u.select_atoms("protein and resid 889 and name CA and segid PROB",updating=True)
res_927_A = u.select_atoms("protein and resid 927 and name CA and segid PROA",updating=True)
res_927_B = u.select_atoms("protein and resid 927 and name CA and segid PROB",updating=True)
# print(L543_A)
df = pd.DataFrame()

with open("TM10_LENGTH_ANGLE_%s.dat"%(systemname),'w+') as ofile:
    ofile.write("#FRAME L_A L_B theta_A theta_B\n")
    for ts in tqdm(u.trajectory[::traj_skip]):
        # d_850 = np.average(ref_850_A.positions,ref_850_B.positions)
        mid_850 = (ref_850_A.positions+ref_850_B.positions)/2
        mid_887 = (ref_887_A.positions+ref_887_B.positions)/2
        unit_z = mid_850-mid_887
        tm10_vec_A = res_889_A.positions-ref_887_A.positions
        tm10_vec_B = res_889_B.positions-ref_887_B.positions
        tm10_length_A = math.dist(*ref_887_A.positions,*res_889_A.positions)
        tm10_length_B = math.dist(*ref_887_B.positions,*res_889_B.positions)
        # print(unit_z[0],tm10_vec_A[0])
    
        angle_A = (np.dot(unit_z[0],tm10_vec_A[0]))/(np.sqrt(np.sum(unit_z[0]**2))*np.sqrt(np.sum(tm10_vec_A[0]**2)))
        angle_A_degree = np.arccos(angle_A)*180/3.14
        angle_B = (np.dot(unit_z[0],tm10_vec_B[0]))/(np.sqrt(np.sum(unit_z[0]**2))*np.sqrt(np.sum(tm10_vec_B[0]**2)))
        angle_B_degree = np.arccos(angle_B)*180/3.14
        ofile.write("%s\t%s\t%s\t%s\t%s\n"%(ts.frame,tm10_length_A,tm10_length_B,angle_A_degree,angle_B_degree))
        ofile.flush()