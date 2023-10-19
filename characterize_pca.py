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
def pca_calculation(selected_str,nComponent,skipping):
    # Aligning the traj to the REF frame. ALIGNMENT ON EACH SEGMENT
    aligner = align.AlignTraj(u,u0,
                          select=selected_str,
                          in_memory=True).run(step=skipping)
    # print(aligner)
    # print("############# Principal Component Analysis\n")
    pc = pca.PCA(u,
                select='protein and backbone',
                # align=True,
                n_components=nComponent).run(step=skipping)
    return pc

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

parser.add_argument('--npc', type=int, default=3,
                    help='Number of Principal Components')

parser.add_argument('--sel', type=str, default='backbone',
                    help='Using VMD selection with " ". Default: backbone')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
traj_begin = args.begin
traj_end = args.end
systemname = args.system
selectionatom=args.sel
ncomponent = args.npc

u = mda.Universe(top_file,traj_file,in_memory=True)
u0 = mda.Universe(top_file,traj_file)
u0.trajectory[0]

n_atom_origin = len(u.atoms)

# Define SELECTIONS
chainA = u.select_atoms("protein and segid PROA and %s"%(selectionatom),updating=True)
chainB = u.select_atoms("protein and segid PROB and %s"%(selectionatom),updating=True)

chain_str = ['protein and segid PROA and %s'%(selectionatom),
             'protein and segid PROB and %s'%(selectionatom)]

chain_list = ['PROA',
              'PROB']

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

with open("PCA_DATA_%s.csv"%(systemname),"w+") as pca_out:
    pca_out.write("# Frame PC1 PC2 PC3 Chain ps\n")
    for ind, (chain) in enumerate((chain_str)):
        pca_result = pca_calculation(chain,ncomponent,skipping=traj_skip)
        selected_segment = u.select_atoms(chain)
        transformed = pca_result.transform(selected_segment, n_components=ncomponent)

        # Projecting PC1 to structure
        print("############# Projecting PC1 to structure\n")
        pc1 = pca_result.p_components[:,0]
        trans1 = transformed[:,0]
        projected = np.outer(trans1, pc1) + pca_result.mean.flatten()
        coordinates = projected.reshape(len(trans1), -1, 3)
        proj1 = mda.Merge(selected_segment)
        proj1.load_new(coordinates, order="fac")
        proj1_selected = proj1.select_atoms("all")
        proj1_selected.write('PCA_PC1_PROJECTED_%s_%s.pdb'%(chain_list[ind],systemname))
        proj1_selected.write('PCA_PC1_PROJECTED_%s_%s.dcd'%(chain_list[ind],systemname),frames='all')

        df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(ncomponent)])
        if ind ==0:
            df['Chain'] = 'A'
        else: df['Chain'] = 'B'
        df['ps'] = df.index * u.trajectory.dt
        for row in tqdm(df.itertuples(),desc='Writing output for %s'%(chain),total=len(u.trajectory)):
            pca_out.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(row[0],row[1],row[2],row[3],row[4],row[5]))