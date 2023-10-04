# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
from MDAnalysis.analysis import hole2
import os 
import warnings
warnings.filterwarnings("ignore")
import re
import mdtraj as md
import sys
# import dask
# import dask.multiprocessing
# dask.config.set(scheduler='processes')
# from dask.distributed import Client
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.visualization import streamlines
import matplotlib.pyplot as plt
import math

top_file = sys.argv[1]   #PSF
traj_file = sys.argv[2]


u = mda.Universe(top_file,traj_file)
# u = u.trajectory[1::100]

boxsize = u.dimensions

chainA = u.select_atoms("segid PROA",updating=True).center_of_geometry()
chainB = u.select_atoms("segid PROB",updating=True).center_of_geometry()

all = u.select_atoms("all")
u.trajectory.next()
all_coord = all.positions
# print(all_coord[:,0])
xmin=all_coord[:,0].min()
xmax=all_coord[:,0].max()
ymin=all_coord[:,1].min()
ymax=all_coord[:,1].max()


# print(xmin,xmax,ymin,ymax)
spacing = 10

u1, v1, average_displacement, standard_deviation_of_displacement = \
    streamlines.generate_streamlines(top_file,
                                    traj_file,
                                    grid_spacing=spacing,
                                    MDA_selection='resname PLPI24 and name P4', 
                                    start_frame=2, end_frame=10000,
                                    xmin=xmin, xmax= xmax,
                                    ymin=ymin, ymax=ymax,
                                    maximum_delta_magnitude=1.0, num_cores=32)

nbinx = u1.shape[1]
nbiny = v1.shape[0]

x = np.linspace(xmin,xmax, nbinx)
y = np.linspace(ymin,ymax, nbiny)

X,Y = np.meshgrid(x,y)
speed = np.sqrt(u1*u1 + v1*v1)

# print(u1.shape)
# print(v1.shape)
# print()
# print(X.shape)
# print(Y.shape)

# print(u1,v1)
# fig = plt.figure()
# ax = fig.add_subplot(111, aspect='equal')
# ax.set_xlabel('x ($\AA$)')
# ax.set_ylabel('y ($\AA$)')
# ax.plot(chainA[0],chainA[1],marker="o")
# ax.text(chainA[0],chainA[1],"CHAIN A")
# ax.plot(chainB[0],chainB[1],marker="o")
# ax.text(chainB[0],chainB[1],"CHAIN B")
# ax.streamplot(X, Y, u1, v1, density=(10,10), color=speed, linewidth=3*speed/speed.max(),cmap='Spectral')

# fig.savefig('streamline.png',dpi=300)

# # # plt.legend()
# plt.show()
