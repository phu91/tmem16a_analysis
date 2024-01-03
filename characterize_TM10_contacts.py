# Updated on 11/29/2022
import MDAnalysis as mda
import numpy as np
import math, sys
import argparse
from MDAnalysis.analysis import align, distances
from MDAnalysis.coordinates.memory import MemoryReader
from tqdm import tqdm
import pandas as pd
def distance_per_residue(selection_tm10,selection_helix):
    selection_tm10_positions = selection_tm10.positions
    selection_tm10_atoms = selection_tm10.atoms

    selection_helix_positions = selection_helix.positions
    selection_helix_atoms = selection_helix.atoms
    
    # selection_helix_positions = u.select_atoms(selection_helix).positions
    # print(selection_tm10_residue,selection_tm10_positions)
    # for i in selection_tm10_residues:
    #     for j,k in zip(selection_tm10_atoms,selection_tm10_positions):
    #         print(i,j,k)
    selection_tm10_label=[]
    for resn,resi in zip(selection_tm10_atoms.resnames,selection_tm10_atoms.resids):
        selection_tm10_label.append((resn+str(resi)))
    # print(len(selection_tm10_label))
    selection_helix_label=[]
    for resn,resi in zip(selection_helix_atoms.resnames,selection_helix_atoms.resids):
        selection_helix_label.append((resn+str(resi)))

    # print(selection_helix_label)
    dist_arr = distances.distance_array(selection_helix_positions,
                                        selection_tm10_positions,
                                    box=u.dimensions)
    # print(len(selection_tm10_atoms),len(selection_helix_atoms),np.shape(dist_arr))
    distance_df = pd.DataFrame(dist_arr,columns=selection_tm10_label)
    distance_df['TM10']=selection_helix_label
    # small.to_csv('tmp.csv')
    distance_df=distance_df.T.reset_index().groupby(['index'],sort=False).min()
    distance_df2 = distance_df.T.groupby(['TM10'],sort=False).min()
    # small2.insert(0, "frame", frame)
    return distance_df2

# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate for Radial Distribution Function calculations. Default = 1 frames')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Add a system name to output file')

parser.add_argument('--cutoff', type=int,default='5',
                    help='Cut-off distance. Default 5A')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
cutoff = args.cutoff
systemname = args.system

u = mda.Universe(top_file,traj_file,in_memory=True)

# Define NTD of TM10

helix_10a = u.select_atoms('segid PROA and resid 887 to 927 and not name H*',updating=True)
helix_10b = u.select_atoms('segid PROB and resid 887 to 927 and not name H*',updating=True)

chaina = u.select_atoms("segid PROA and not name H*",updating=True)
chainb = u.select_atoms("segid PROB and not name H*",updating=True)

# print(helix_10a.atoms)

with open ("DISTANCE_TM10_PROFILE_%s"%(systemname),'w+') as ofile:
    for ts in u.trajectory[::traj_skip]:
        distance_profile_a = distance_per_residue(helix_10a,chainb)
        distance_profile_a.insert(0,'FRAME',ts.frame)
        print("FRAME %s"%(ts.frame))
        # print(distance_profile_a)
        distance_profile_a.to_csv("DISTANCE_TM10_PROFILE_%s"%(systemname),mode="a",index=True,header=True)

        # distance_profile_b = distance_per_residue(helix_10b,chaina)
        # # print(len(distance_profile_a),len(distance_profile_b))
        # distance_profile_a=np.array(distance_profile_a)
        # distance_profile_b=np.array(distance_profile_b)
        # for x in tqdm(distance_profile_a,desc='Writing chain A Frame %s'%(ts.frame)):
        #     ofile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(ts.frame,x[0],x[1],x[2],x[3],x[4],x[5],x[6]))        
        #     ofile.flush()
        # for x in tqdm(distance_profile_b,desc='Writing chain B Frame %s'%(ts.frame)):
        #     ofile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(ts.frame,x[0],x[1],x[2],x[3],x[4],x[5],x[6]))        
        #     ofile.flush()