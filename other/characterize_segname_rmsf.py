# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
import math, sys
import argparse
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import rms, align
import warnings
from tqdm import tqdm
warnings.filterwarnings('ignore')

# FUNCTIONS
def rmsf_calculation(selected_str,skipping):
    # ref = u.select_atoms(selected_str)      ## FIRST FRAME
    # Aligning the traj to the REF frame. ALIGNMENT ON EACH SEGMENT
    # print(len(u.atoms))
    average = align.AverageStructure(u, u, select=selected_str,
                                 ref_frame=0).run(step=skipping)
    print("Aligning TRAJ: %s  || to the AVERAGED structure || %s"%(selected_str,selected_str))
    ref = average.results.universe
    aligner = align.AlignTraj(u, ref,
                          select=selected_str,
                          in_memory=True).run(step=skipping)

    selected_segment = u.select_atoms(selected_str)
    R = rms.RMSF(selected_segment, verbose=True).run(step=skipping)
    return R

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

u = mda.Universe(top_file,traj_file,in_memory=True)
u.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
n_atom_origin = len(u.atoms)

# Define REFERENCE: FRAME 0
chainA = u.select_atoms("protein and segid PROA and name CA",updating=True)
chainB = u.select_atoms("protein and segid PROB and name CA",updating=True)
chainI = u.select_atoms("protein and segid PROI and name CA",updating=True)

chainA_full = u.select_atoms("protein and segid PROA",updating=True)
chainB_full = u.select_atoms("protein and segid PROB",updating=True)
chainI_full = u.select_atoms("protein and segid PROI",updating=True)

# print(chainA.residues)
chain_str = ['protein and segid PROA and name CA','protein and segid PROB and name CA','protein and segid PROI and name CA']
ref_list = [chainA_full,chainB_full,chainI_full]
chain_list = ['LTGFb1A','LTGFb1B','GARP']

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
    print("########################################################\n")
else:
    print("\n########################################################")
    print("TOP : %s"%(top_file))
    print("TRAJ: %s"%(traj_file))
    print("NUMBER OF ATOMS (ORIGINAL): %s"%(n_atom_origin))
    print("NUMBER OF ATOMS (CURRENT) : %s"%(len(u.atoms)))
    print("########################################################\n")
    pass

with open("RMSF_%s.dat"%(systemname),"w+") as rmsf_out:
    for ind, (chain,ref) in enumerate(zip(chain_str,ref_list)):
        rmsf_out.write("#chain resname resid rmsf sysname\n")
        RMSF = rmsf_calculation(chain,traj_skip)
        for RES,rmsf in tqdm(zip(ref.residues,RMSF.results.rmsf),total=len(ref.residues),desc=chain):
            # print(RES,chain,rmsf)
            if ind ==0:
                # print("A",RES.resid)
                rmsf_out.write("%s\t%s\t%s\t%s\t%s\n"%(chain_list[ind],RES.resname,RES.resid,rmsf,systemname))
                rmsf_out.flush()
                RES.atoms.tempfactors = rmsf
                ref.atoms.write('rmsf_tempfactors_chain_%s.pdb'%(chain_list[ind]))
            if ind ==1:
                # print("B",RES.resid)
                rmsf_out.write("%s\t%s\t%s\t%s\t%s\n"%(chain_list[ind],RES.resname,RES.resid,rmsf,systemname))
                rmsf_out.flush()
                RES.atoms.tempfactors = rmsf
                ref.atoms.write('rmsf_tempfactors_chain_%s.pdb'%(chain_list[ind]))
            if ind ==2:
                # print("I",RES.resid)
                rmsf_out.write("%s\t%s\t%s\t%s\t%s\n"%(chain_list[ind],RES.resname,RES.resid,rmsf,systemname))
                rmsf_out.flush()
                RES.atoms.tempfactors = rmsf
                ref.atoms.write('rmsf_tempfactors_chain_%s.pdb'%(chain_list[ind]))
