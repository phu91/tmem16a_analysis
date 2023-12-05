# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
import math, sys
import argparse
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import pca, align
from tqdm import tqdm

import warnings
# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

# FUNCTIONS


def extract_frames(startFrame,endFrame):
    whole_backbone = u.select_atoms("protein and backbone")
    with mda.Writer("SEGMENT_FRAME_%s_TO_%s.pdb"%(startFrame,endFrame)) as pdb:
        pdb.write(whole_backbone)

    with mda.Writer("SEGMENT_FRAME_%s_TO_%s.xtc"%(startFrame,endFrame), all.n_atoms) as W:
        for ts in u.trajectory[startFrame:endFrame:traj_skip]:
            W.write(whole_backbone)

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

parser.add_argument('--sel', type=str, default='name CA',
                    help='Using VMD selection with " ". Default: backbone')

parser.add_argument('--ref', type=str, default='step5_input.pdb',
                    help='If not provided, use default output (step5_input.pdb) from CHARMM-GUI as the reference (Default)')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
traj_begin = args.begin
traj_end = args.end
systemname = args.system
selectionatom=args.sel
reference = args.ref

u = mda.Universe(top_file,traj_file,in_memory=True)
u0 = mda.Universe(top_file,reference)


n_atom_origin = len(u.atoms)

# Define SELECTIONS ## USE CA atom for faster testing. Change back to 'backbone' in production run
chainA = u.select_atoms("protein and segid PROA and %s"%(selectionatom),updating=True)
chainB = u.select_atoms("protein and segid PROB and %s"%(selectionatom),updating=True)

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

# with open("TILT_%s.csv"%(systemname),"w+") as pca_out:
#     pca_out.write("# Frame rad\n")

