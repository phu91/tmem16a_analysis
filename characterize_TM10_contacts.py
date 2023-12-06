# Updated on 11/29/2022
import MDAnalysis as mda
import numpy as np
import math, sys
import argparse
from MDAnalysis.analysis import align, distances
from MDAnalysis.coordinates.memory import MemoryReader
from tqdm import tqdm

def distance_per_residue(selection_1,selection_2):
    selection_1_com = selection_1.center_of_mass(compound='residues')
    selection_2_com = selection_2.center_of_mass(compound='residues')
    # print(frame,selection_1_com)
    # print(u.dimensions)
    dist_arr = distances.distance_array(selection_1_com, # reference
                                    selection_2_com, # configuration
                                    box=u.dimensions)
    distance_collect=[]
    for ix,x in enumerate(selection_1.residues):
        for iy,y in enumerate(selection_2.residues):
            # if dist_arr[ix,iy]<=cutoff:
            if x.segid=='PROA':
                distance_collect.append((x.resname,x.resid,'A',y.resname,y.resid,'B',dist_arr[ix,iy]))
            else:
                distance_collect.append((x.resname,x.resid,'B',y.resname,y.resid,'A',dist_arr[ix,iy]))
            # else:
            #     distance_collect.append(('0'))
    return distance_collect
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

helix_10a_str = 'segid PROA and resid 887 to 927'
helix_10b_str = 'segid PROB and resid 887 to 927'

helix_10a = u.select_atoms(helix_10a_str,updating=True)
helix_10b = u.select_atoms(helix_10b_str,updating=True)

chaina = u.select_atoms("segid PROA",updating=True)
chainb = u.select_atoms("segid PROB",updating=True)

with open ("DISTANCE_TM10_PROFILE_%s"%(systemname),'w+') as ofile:
    ofile.write("#frame resname1 resid1 chain1 resname2 resid2 chain2 distance\n")
    for ts in u.trajectory[::traj_skip]:
        distance_profile_a = distance_per_residue(helix_10a,chainb)
        distance_profile_b = distance_per_residue(helix_10b,chaina)
        # print(len(distance_profile_a),len(distance_profile_b))
        distance_profile_a=np.array(distance_profile_a)
        distance_profile_b=np.array(distance_profile_b)
        for x in tqdm(distance_profile_a,desc='Writing chain A Frame %s'%(ts.frame)):
            ofile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(ts.frame,x[0],x[1],x[2],x[3],x[4],x[5],x[6]))        
            ofile.flush()
        for x in tqdm(distance_profile_b,desc='Writing chain B Frame %s'%(ts.frame)):
            ofile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(ts.frame,x[0],x[1],x[2],x[3],x[4],x[5],x[6]))        
            ofile.flush()