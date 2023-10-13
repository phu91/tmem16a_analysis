# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.visualization import streamlines
import matplotlib.pyplot as plt
import math, sys
from MDAnalysis.lib.distances import distance_array
import seaborn as sns 
import argparse

# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate for Radial Distribution Function calculations. Default = 1 frames')

args = parser.parse_args()

top_file =  args.top
traj_file = args.traj
traj_skip = args.skip

u = mda.Universe(top_file,traj_file)

chainA_com = u.select_atoms("segid PROA",updating=True).center_of_geometry()
chainB_com = u.select_atoms("segid PROB",updating=True).center_of_geometry()

pip2 = u.select_atoms("resname PLPI* and name P",updating=True)

# print(chainA_com,chainB_com,pip2)
with open("pip2_track_all.dat",'w+') as of:
    of.write("#frame pip2_id pip2_x pip2_y\n")
    for ts in u.trajectory[::traj_skip]:
        print(ts.frame)
        # distance_a =distance_array(reference=chainA_com, 
        #                            configuration=pip2.positions, 
        #                            box=u.dimensions, 
        #                            result=None, 
        #                            backend='OpenMP'
        # )
        id_pip2 = pip2.resids
        x_pip2 = pip2.positions[:,0]
        y_pip2 = pip2.positions[:,1]
        for i in range(len(x_pip2)):
            of.write("%s\t%s\t%s\t%s\n"%(ts.frame,id_pip2[i],x_pip2[i],y_pip2[i]))
            of.flush()
