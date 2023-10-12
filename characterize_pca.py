# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math, sys
import seaborn as sns 
import argparse
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import pca, align
import nglview as nv
from tqdm import tqdm
from nglview.contrib.movie import MovieMaker

import warnings
# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

sns.set_context("paper")

# FUNCTIONS
def pca_calculation(selected_str,nComponent,skipping):
    # Aligning the traj to the REF frame. ALIGNMENT ON EACH SEGMENT
    # aligner = align.AlignTraj(u, u,
    #                       select=selected_str,
    #                       in_memory=True).run(step=skipping)

    pc = pca.PCA(u, 
                select=selected_str,
                align=True,
                n_components=nComponent).run(step=skipping)
    return pc

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

parser.add_argument('--npc', type=int, default=3,
                    help='Number of Principal Components')


args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
traj_begin = args.begin
traj_end = args.end
systemname = args.system
ncomponent = args.npc

u = mda.Universe(top_file,traj_file,in_memory=True)
n_atom_origin = len(u.atoms)

# Define REFERENCE: FRAME 0
chainA = u.select_atoms("protein and segid PROA and backbone",updating=True)
chainB = u.select_atoms("protein and segid PROB and backbone",updating=True)

# chainA_full = u.select_atoms("protein and segid PROA and backbone",updating=True)
# chainB_full = u.select_atoms("protein and segid PROB and backbone",updating=True)

chain_str = ['protein and segid PROA and backbone','protein and segid PROB and backbone']

chain_list = ['PROA','PROB']

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

n_frames = len(u.trajectory)

with open("PCA_TAIL_%s.csv"%(systemname),"w+") as pca_out:
    pca_out.write("# Frame PC1 PC2 PC3 Chain ps\n")
    for ind, (chain) in enumerate((chain_str)):
        pca_result = pca_calculation(chain,ncomponent,skipping=traj_skip)
        selected_segment = u.select_atoms(chain)
        transformed = pca_result.transform(selected_segment, n_components=ncomponent)
        # Projecting PC1 to structure
        pc1 = pca_result.p_components[:,0]
        trans1 = transformed[:,0]
        projected = np.outer(trans1, pc1) + pca_result.mean.flatten()
        coordinates = projected.reshape(len(trans1), -1, 3)
        proj1 = mda.Merge(selected_segment)
        proj1.load_new(coordinates, order="fac")
        proj1_selected = proj1.select_atoms("all")
        proj1_selected.write('TEST_chain_%s.pdb'%(chain_list[ind]))
        proj1_selected.write('TEST_chain_%s.dcd'%(chain_list[ind]),frames='all')
        # print(proj1.trajectory)
        # view = nv.show_mdanalysis(proj1.atoms)
        # print(view)
        # movie = MovieMaker(view,
        #                     step=1,  # keep every 4th frame
        #                     output='pc1.',
        #                     render_params={"factor": 3},  # set to 4 for highest quality
        #                     )
        # proj1.write("c-alpha.dcd")

        df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(ncomponent)])
        if ind ==0: 
            df['Chain'] = 'A'
        else: df['Chain'] = 'B'
        df['ps'] = df.index * u.trajectory.dt
        for row in tqdm(df.itertuples(),desc='Writing output for %s'%(chain),total=n_frames):
            pca_out.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(row[0],row[1],row[2],row[3],row[4],row[5]))