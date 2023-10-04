# Updated on 11/29/2022
import pandas as pd
import numpy as np
import math, sys
import argparse
import mdtraj as md
import MDAnalysis as mda
from multiprocessing import Pool, Process,Value, Array
# FUNCTIONS

def dssp_calculation(frame):
    with open("DSSP_%s.dat"%(top_file.rstrip('.pdbsf')),"w+") as dssp_out:
        dssp_out.write("#frame helix chain helicity\n")
        H_content = []
        for ind, (helixa,helixb) in enumerate(zip(helixa_list,helixb_list)):
            current_frame = md.load_frame(traj_file,frame,top=top_file,atom_indices=helixa.indices)
            dssp1 = md.compute_dssp(current_frame)
            dssp1 = dssp1.tolist()
            for i in range(len(dssp1)):
                n_residues = len(dssp1[i])
                count_H = dssp1[i].count('H')
                percent_H = count_H/n_residues
                # return (frame,ind+1,"A",percent_H)
                # print(ind)
                # H_content = np.append(H_content,(frame,ind+1,"A",percent_H))
                # dssp_out.write("%s\t%s\t%s\t%s\n"%(frame,ind+1,"A",percent_H))
        return (H_content)
        # H = H_content.tolist()
        # return (H)
        # dssp_out.write("%s\t%s\t%s\t%s\n"%(H_content[0],H_content[1],H_content[2],H_content[3]))
        # dssp_out.flush()
        # print(*H_content)
        # return H_content

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

# parser.add_argument('--chain', type=str,
#                     help='Chain of protein')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
systemname = args.system
# chain_name = args.chain


u = mda.Universe(top_file,traj_file)

# Define HELIX
helix_01a = u.select_atoms("segid PROA and resid 327 to 360")
helix_02a = u.select_atoms("segid PROA and resid 399 to 439")
helix_03a = u.select_atoms("segid PROA and resid 493 to 518")
helix_04a = u.select_atoms("segid PROA and resid 534 to 562")
helix_05a = u.select_atoms("segid PROA and resid 572 to 598")
helix_06a = u.select_atoms("segid PROA and resid 629 to 664")
helix_07a = u.select_atoms("segid PROA and resid 692 to 713")
helix_08a = u.select_atoms("segid PROA and resid 714 to 738")
helix_09a = u.select_atoms("segid PROA and resid 753 to 775")
helix_10a = u.select_atoms("segid PROA and resid 854 to 881")

helix_01b = u.select_atoms("segid PROB and resid 327 to 360")
helix_02b = u.select_atoms("segid PROB and resid 399 to 439")
helix_03b = u.select_atoms("segid PROB and resid 493 to 518")
helix_04b = u.select_atoms("segid PROB and resid 534 to 562")
helix_05b = u.select_atoms("segid PROB and resid 572 to 598")
helix_06b = u.select_atoms("segid PROB and resid 629 to 664")
helix_07b = u.select_atoms("segid PROB and resid 692 to 713")
helix_08b = u.select_atoms("segid PROB and resid 714 to 738")
helix_09b = u.select_atoms("segid PROB and resid 753 to 775")
helix_10b = u.select_atoms("segid PROB and resid 854 to 881")

helixa_list = [helix_01a,
helix_02a,
helix_03a,
helix_04a,
helix_05a,
helix_06a,
helix_07a,
helix_08a,
helix_09a,
helix_10a]

helixb_list = [helix_01b,
helix_02b,
helix_03b,
helix_04b,
helix_05b,
helix_06b,
helix_07b,
helix_08b,
helix_09b,
helix_10b]

helixa_list_label = ['helix_01a',
'helix_02a',
'helix_03a',
'helix_04a',
'helix_05a',
'helix_06a',
'helix_07a',
'helix_08a',
'helix_09a',
'helix_10a']

helixb_list_label = ['helix_01b',
'helix_02b',
'helix_03b',
'helix_04b',
'helix_05b',
'helix_06b',
'helix_07b',
'helix_08b',
'helix_09b',
'helix_10b']

frame_list = np.arange(0,u.trajectory.n_frames,traj_skip)
#     dssp_out.write("#frame helix chain helicity\n")

if __name__ == '__main__':
    with Pool(10) as p:
        results = p.map_async(dssp_calculation, frame_list)
        # results.wait()
        p.close()
        p.join()
        print(results.get())

    # for result in results:
        # print(result)


# with open("DSSP_%s.dat"%(top_file.rstrip('.pdbsf')),"w+") as dssp_out:
#     dssp_out.write("#frame helix chain helicity\n")
#     for ts in u.trajectory[::traj_skip]:
#         for ind, (helixa,helixb) in enumerate(zip(helixa_list,helixb_list)):
#            Ha = dssp_calculation(top_file,traj_file,ts.frame,helixa.indices)
#            Hb = dssp_calculation(top_file,traj_file,ts.frame,helixb.indices)
#            dssp_out.write("%s\t%s\t%s\t%s\n"%(ts.frame,ind+1,"A",*Ha))
#            dssp_out.write("%s\t%s\t%s\t%s\n"%(ts.frame,ind+1,"B",*Hb))
#            dssp_out.flush()