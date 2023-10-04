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
def rmsf_ref(selected_str,skipping):
    average = align.AverageStructure(u, u, select=selected_str,
                                 ref_frame=0).run(step=skipping)
    ref = average.results.universe
    # Aligning the traj to the REF frame. ALIGNMENT ON EACH SEGMENT
    print("Aligning TRAJ: %s to the AVERAGED structure: %s\n"%(selected_str,selected_str))
    aligner = align.AlignTraj(u, ref,
                          select=selected_str,
                          in_memory=True).run(step=skipping)

    selected_segment = u.select_atoms(selected_str)
    R = RMSF(selected_segment, verbose=True).run(step=skipping)
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

parser.add_argument('--systemid', type=str,
                    help='ID of the SYSTEM')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
systemname = args.system
systemid = args.systemid


u = mda.Universe(top_file,traj_file,in_memory=True)

# Define SEGMENT
chainA = u.select_atoms("protein and segid PROA and name CA",updating=True)
chainB = u.select_atoms("protein and segid PROB and name CA",updating=True)

print(chainA.residues)
chain_str = ['protein and segid PROA and name CA','protein and segid PROB and name CA']
chain_list = [chainA,chainB]

with open("RMSF_SEGid_%s.dat"%(systemname),"w+") as rmsf_out:
    rmsf_out.write("#chain resid rmsf sys\n")
#     # Average Structure
    for ind, (chain,sel) in enumerate(zip(chain_str,chain_list)):
        rmsf = rmsf_ref(chain,traj_skip)

        u.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
        for residue, r_value in zip(sel.residues, rmsf.results.rmsf):
            residue.atoms.tempfactors = r_value
        if ind ==0:
            u.atoms.write('rmsf_segidA_tempfactors.pdb')
        else:
            u.atoms.write('rmsf_segidA_tempfactors.pdb')

        for resid,rmsf in zip(sel.resids,rmsf.results.rmsf):
            if ind ==0:
                print("A",resid,rmsf)
                rmsf_out.write("%s\t%s\t%s\t%s\n"%("A",resid,rmsf,systemid))
                rmsf_out.flush()
            else:
                print("B",resid,rmsf)
                rmsf_out.write("%s\t%s\t%s\t%s\n"%("B",resid,rmsf,systemid))
                rmsf_out.flush()

