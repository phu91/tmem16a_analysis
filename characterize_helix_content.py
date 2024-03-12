# Updated on 11/29/2022
import pandas as pd
import numpy as np
import math, sys
import argparse
import mdtraj as md
import MDAnalysis as mda
from multiprocessing import Pool, Process,Value, Array
# FUNCTIONS

# def dssp_calculation(frame):
#     with open("DSSP_%s.dat"%(top_file.rstrip('.pdbsf')),"w+") as dssp_out:
#         dssp_out.write("#frame helix chain helicity\n")
#         H_content = []
#         for ind, (helixa,helixb) in enumerate(zip(helixa_list,helixb_list)):
#             current_frame = md.load_frame(traj_file,frame,top=top_file,atom_indices=helixa.indices)
#             dssp1 = md.compute_dssp(current_frame)
#             dssp1 = dssp1.tolist()
#             for i in range(len(dssp1)):
#                 n_residues = len(dssp1[i])
#                 count_H = dssp1[i].count('H')
#                 percent_H = count_H/n_residues
#                 # return (frame,ind+1,"A",percent_H)
#                 # print(ind)
#                 # H_content = np.append(H_content,(frame,ind+1,"A",percent_H))
#                 # dssp_out.write("%s\t%s\t%s\t%s\n"%(frame,ind+1,"A",percent_H))
#         return (H_content)
        # H = H_content.tolist()
        # return (H)
        # dssp_out.write("%s\t%s\t%s\t%s\n"%(H_content[0],H_content[1],H_content[2],H_content[3]))
        # dssp_out.flush()
        # print(*H_content)
        # return H_content

###############
def extract_frames(startFrame,endFrame,selections):
    protein = u.select_atoms(selections)
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

parser.add_argument('--sel', type=str, default='all',
                    help='Add a selections of atoms to extract. Default = "all"')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
traj_begin = args.begin
traj_end = args.end
systemname = args.system
selected_extraction=args.sel

u = mda.Universe(top_file,traj_file)   #,in_memory=True)
n_atom_origin = len(u.atoms)
print(u)
### SPLICING TRAJECTORY
if traj_end != -1:
    extract_frames(traj_begin,traj_end,selected_extraction)
    u = mda.Universe("SEGMENT_FRAME_%s_TO_%s.pdb"%(traj_begin,traj_end),
                     "SEGMENT_FRAME_%s_TO_%s.xtc"%(traj_begin,traj_end),
                     in_memory=True)
    print("\n########################################################")
    print("TOP : SEGMENT_FRAME_%s_TO_%s.pdb"%(traj_begin,traj_end))
    print("TRAJ: SEGMENT_FRAME_%s_TO_%s.xtc"%(traj_begin,traj_end))
    print("NUMBER OF ATOMS (ORIGINAL): %s"%(n_atom_origin))
    print("NUMBER OF ATOMS (CURRENT) : %s"%(len(u.atoms)))
    print("NUMBER OF FRAME (CURRENT) : %s"%(len(u.trajectory)))

    print("########################################################\n")
else:
    print("\n########################################################")
    print("TOP : %s"%(top_file))
    print("TRAJ: %s"%(traj_file))
    print("NUMBER OF ATOMS (ORIGINAL): %s"%(n_atom_origin))
    print("NUMBER OF ATOMS (CURRENT) : %s"%(len(u.atoms)))
    print("NUMBER OF FRAME (CURRENT) : %s"%(len(u.trajectory)))

    print("########################################################\n")
    pass


### Define HELIX UPDATED BASED ON BENDIX
# helix_01a = u.select_atoms("segid PROA and resid 327 to 360",updating=True)
# helix_02a = u.select_atoms("segid PROA and resid 399 to 439",updating=True)
# helix_03a = u.select_atoms("segid PROA and resid 478 to 518",updating=True)
# helix_04a = u.select_atoms("segid PROA and resid 534 to 562",updating=True)
# helix_05a = u.select_atoms("segid PROA and resid 572 to 598",updating=True)
helix_06a = u.select_atoms("segid PROA and resid 637 to 648",updating=True)
# helix_07a = u.select_atoms("segid PROA and resid 692 to 713",updating=True)
# helix_08a = u.select_atoms("segid PROA and resid 718 to 740",updating=True)
# helix_09a = u.select_atoms("segid PROA and resid 753 to 780",updating=True)
# helix_10a = u.select_atoms("segid PROA and resid 887 to 925",updating=True)

# helix_01b = u.select_atoms("segid PROB and resid 327 to 360",updating=True)
# helix_02b = u.select_atoms("segid PROB and resid 399 to 439",updating=True)
# helix_03b = u.select_atoms("segid PROB and resid 478 to 518",updating=True)
# helix_04b = u.select_atoms("segid PROB and resid 534 to 562",updating=True)
# helix_05b = u.select_atoms("segid PROB and resid 572 to 598",updating=True)
helix_06b = u.select_atoms("segid PROB and resid 634 to 648",updating=True)
# helix_07b = u.select_atoms("segid PROB and resid 692 to 713",updating=True)
# helix_08b = u.select_atoms("segid PROB and resid 718 to 740",updating=True)
# helix_09b = u.select_atoms("segid PROB and resid 753 to 780",updating=True)
# helix_10b = u.select_atoms("segid PROB and resid 887 to 925",updating=True)


helixa_list = [
# helix_01a,
# helix_02a,
# helix_03a,
# helix_04a,
# helix_05a,
helix_06a,
# helix_07a,
# helix_08a,
# helix_09a,
# helix_10a
]

helixb_list = [
# helix_01b,
# helix_02b,
# helix_03b,
# helix_04b,
# helix_05b,
helix_06b,
# helix_07b,
# helix_08b,
# helix_09b,
# helix_10b
]

# helixa_list_label = ['helix_01a',
# 'helix_02a',
# 'helix_03a',
# 'helix_04a',
# 'helix_05a',
# 'helix_06a',
# 'helix_07a',
# 'helix_08a',
# 'helix_09a',
# 'helix_10a']

# helixb_list_label = ['helix_01b',
# 'helix_02b',
# 'helix_03b',
# 'helix_04b',
# 'helix_05b',
# 'helix_06b',
# 'helix_07b',
# 'helix_08b',
# 'helix_09b',
# 'helix_10b']

# print(u)
frame_list = np.arange(0,u.trajectory.n_frames,traj_skip)
#     dssp_out.write("#frame helix chain helicity\n")

with open("DSSP_%s.dat"%(top_file.rstrip('.pdbsf')),"w+") as dssp_out:
    dssp_out.write("#frame resname resid helix_a helix_b\n")
    for ts in u.trajectory[traj_begin:traj_end:traj_skip]:
        for ind, (helixa,helixb) in enumerate(zip(helixa_list,helixb_list)):
            # print(len(helixa.residues))   ## 0-based index used in mdtraj and MDAnalysis
            current_frame_a = md.load_frame(traj_file,ts.frame,top=top_file,atom_indices=helixa.indices)
            dssp_a = md.compute_dssp(current_frame_a,simplified=False)
            current_frame_b = md.load_frame(traj_file,ts.frame,top=top_file,atom_indices=helixb.indices)
            dssp_b = md.compute_dssp(current_frame_b,simplified=False)

            for i,j,k in zip(helixa.residues,*dssp_a,*dssp_b):
                # print(j,len(j))
                dssp_out.write("%s\t%s\t%s\t%s\t%s\n"%(ts.frame,i.resname,i.resid,j,k))
                dssp_out.flush()
        print(ts.frame)
