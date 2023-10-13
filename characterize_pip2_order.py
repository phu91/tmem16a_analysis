# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.visualization import streamlines
import math, sys
from MDAnalysis.lib.distances import distance_array
import argparse
from tqdm import tqdm

### FUNCTION
def lipid_enrichment(helix):
    count=0
    distances_a=distance_array(reference=pip2.positions, 
                           configuration=helix.positions, 
                           box=u.dimensions, 
                           result=None, 
                           backend='OpenMP'
    )
    for resOnHelix in distances_a:
        for distance in resOnHelix:
            if distance <=rcutoff:
                # print(resOnHelix)
                count+=1
                break    ## break the loop when one distance is less than rcut

    # print(count)
    return count




# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

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

parser.add_argument('--rcut', type=int, default=3.5,
                    help='R cutoff value. Default = 3.5 Angstrom')

args = parser.parse_args()

top_file =  args.top
traj_file = args.traj
traj_begin = args.begin
traj_end = args.end
traj_skip = args.skip
systemname = args.system
rcutoff = args.rcut

u = mda.Universe(top_file,traj_file)

chainA_com = u.select_atoms("segid PROA",updating=True)
chainB_com = u.select_atoms("segid PROB",updating=True)

sn1 = u.select_atoms("name C2 C21 C210 C211 C212 C213 C214 C215 C216 C217 C218 C22 C23 C24 C25 C26 C27 C28 C29",updating=True)
sn2 = u.select_atoms("name C3 C31 C310 C311 C312 C313 C314 C315 C316 C317 C318 C32 C33 C34 C35 C36 C37 C38 C39",updating=True)

print(len(sn1))

# print("\n###### R_cutoff is {} Angstrom\n".format(rcutoff))

# print(chainA_com,chainB_com,pip2)
# with open("PIP2_ORDERED_%s.dat"%(systemname),'w+') as of:
#     of.write("#frame chain helix npip2\n")
#     for ts in tqdm(u.trajectory[traj_begin:traj_end:traj_skip]):
#         for ind, (helixa,helixb) in enumerate(zip(helixa_list,helixb_list)):
#             counta = lipid_enrichment(helixa)
#             countb = lipid_enrichment(helixb)

# scc = SCC(
#   universe=u,
#   tail_sel="name ??A"
# )