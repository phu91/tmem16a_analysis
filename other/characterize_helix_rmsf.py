# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math, sys
import seaborn as sns 
import argparse
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.coordinates.memory import MemoryReader

sns.set_context("paper")

# FUNCTIONS
def alignment_ref(helix_str,skipping):
    average = align.AverageStructure(u, u, select=helix_str,
                                 ref_frame=0).run(step=skipping)
    ref = average.results.universe
    # Aligning the traj to the REF frame
    print("Aligning %s to the REF\n"%(helix_str))
    aligner = align.AlignTraj(u, ref,
                          select=helix_str,
                          in_memory=True).run(step=skipping)

def rmsf_helix_analysis(selected_helix_str,skipping):
    selected_helix = u.select_atoms(selected_helix_str+" and name CA")
    R = RMSF(selected_helix, verbose=True).run(step=skipping)
    return R

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


u = mda.Universe(top_file,traj_file,in_memory=True)

# Define HELIX
helix_01a_str = 'segid PROA and resid 327 to 360 and name CA'
helix_02a_str = 'segid PROA and resid 399 to 439 and name CA'
helix_03a_str = 'segid PROA and resid 493 to 518 and name CA'
helix_04a_str = 'segid PROA and resid 534 to 562 and name CA'
helix_05a_str = 'segid PROA and resid 572 to 598 and name CA'
helix_06a_str = 'segid PROA and resid 629 to 664 and name CA'
helix_07a_str = 'segid PROA and resid 692 to 713 and name CA'
helix_08a_str = 'segid PROA and resid 714 to 738 and name CA'
helix_09a_str = 'segid PROA and resid 753 to 775 and name CA'
helix_10a_str = 'segid PROA and resid 854 to 881 and name CA'

helix_01b_str = 'segid PROB and resid 327 to 360 and name CA'
helix_02b_str = 'segid PROB and resid 399 to 439 and name CA'
helix_03b_str = 'segid PROB and resid 493 to 518 and name CA'
helix_04b_str = 'segid PROB and resid 534 to 562 and name CA'
helix_05b_str = 'segid PROB and resid 572 to 598 and name CA'
helix_06b_str = 'segid PROB and resid 629 to 664 and name CA'
helix_07b_str = 'segid PROB and resid 692 to 713 and name CA'
helix_08b_str = 'segid PROB and resid 714 to 738 and name CA'
helix_09b_str = 'segid PROB and resid 753 to 775 and name CA'
helix_10b_str = 'segid PROB and resid 854 to 881 and name CA'

helix_01a = u.select_atoms("segid PROA and resid 327 to 360 and name CA")
helix_02a = u.select_atoms("segid PROA and resid 399 to 439 and name CA")
helix_03a = u.select_atoms("segid PROA and resid 493 to 518 and name CA")
helix_04a = u.select_atoms("segid PROA and resid 534 to 562 and name CA")
helix_05a = u.select_atoms("segid PROA and resid 572 to 598 and name CA")
helix_06a = u.select_atoms("segid PROA and resid 629 to 664 and name CA")
helix_07a = u.select_atoms("segid PROA and resid 692 to 713 and name CA")
helix_08a = u.select_atoms("segid PROA and resid 714 to 738 and name CA")
helix_09a = u.select_atoms("segid PROA and resid 753 to 775 and name CA")
helix_10a = u.select_atoms("segid PROA and resid 854 to 881 and name CA")
helix_01b = u.select_atoms("segid PROB and resid 327 to 360 and name CA")
helix_02b = u.select_atoms("segid PROB and resid 399 to 439 and name CA")
helix_03b = u.select_atoms("segid PROB and resid 493 to 518 and name CA")
helix_04b = u.select_atoms("segid PROB and resid 534 to 562 and name CA")
helix_05b = u.select_atoms("segid PROB and resid 572 to 598 and name CA")
helix_06b = u.select_atoms("segid PROB and resid 629 to 664 and name CA")
helix_07b = u.select_atoms("segid PROB and resid 692 to 713 and name CA")
helix_08b = u.select_atoms("segid PROB and resid 714 to 738 and name CA")
helix_09b = u.select_atoms("segid PROB and resid 753 to 775 and name CA")
helix_10b = u.select_atoms("segid PROB and resid 854 to 881 and name CA")

helixa_str_list = [helix_01a_str,
helix_02a_str,
helix_03a_str,
helix_04a_str,
helix_05a_str,
helix_06a_str,
helix_07a_str,
helix_08a_str,
helix_09a_str,
helix_10a_str]

helixb_str_list = [helix_01b_str,
helix_02b_str,
helix_03b_str,
helix_04b_str,
helix_05b_str,
helix_06b_str,
helix_07b_str,
helix_08b_str,
helix_09b_str,
helix_10b_str]

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

with open("RMSF_profile_%s.dat"%(systemname),"w+") as rdf_out:
    rdf_out.write("#helix chain resid rmsf\n")
    # Average Structure
    for ind, (helixa,helixb,helixa_str,helixb_str) in enumerate(zip(helixa_list,helixb_list,helixa_str_list,helixb_str_list)):
        alignment_ref(helixa_str,traj_skip)
        Ra = rmsf_helix_analysis(helixa_str,traj_skip)

        alignment_ref(helixb_str,traj_skip)
        Rb = rmsf_helix_analysis(helixb_str,traj_skip)
        for resid,rmsf in zip(helixa.resids,Ra.results.rmsf):
            print(ind+1,"A",resid,rmsf)
            rdf_out.write("%s\t%s\t%s\t%s\n"%(ind+1,"A",resid,rmsf))
            rdf_out.flush()
        for resid,rmsf in zip(helixb.resids,Rb.results.rmsf):
            print(ind+1,"B",resid,rmsf)
            rdf_out.write("%s\t%s\t%s\t%s\n"%(ind+1,"B",resid,rmsf))
            rdf_out.flush()


# def dssp_calculation(top_file,traj_file,frame):
#     current_frame = md.load(traj_file,top=top_file)
#     # print(current_frame)
#     dssp1 = md.compute_dssp(current_frame)
#     dssp1 = dssp1.tolist()
#     H_content = []
#     for i in range(len(dssp1)):
#         n_residues = len(dssp1[i])
#         count_H = dssp1[i].count('H')
#         percent_H = count_H/n_residues
#         H_content = np.append(H_content,percent_H)
    # return H_content

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