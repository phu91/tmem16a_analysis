# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
import os 
import warnings
warnings.filterwarnings("ignore")
import re
import mdtraj as md
import sys
import argparse

from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import hole2


def pore_size(pore_center,frame,chain):
    # try:
    # 	os.mkdir("./hole2_output_%s"%(chain))
    # except OSError as error:
    # 	print(error)
    os.makedirs("./hole2_output_%s"%(chain),exist_ok=True)
    working_dir ='hole2_output_%s'%(chain)
    sel = u.select_atoms("protein")

    sel.write(working_dir+"/tmp.pdb") 
    profiles = hole2.hole(working_dir+"/tmp.pdb",
                    executable='~/hole2/exe/hole',
                    outfile=working_dir+'/hole%i.out'%(frame),
                    sphpdb_file=working_dir+'/hole%i.sph'%(frame),
                    vdwradii_file=None,
                    random_seed=31415,
                    keep_files=True,
                    cpoint=pore_center,
                    cvect=[0,0,1])

    # hole2.create_vmd_surface(filename='hole2_output/hole%i.vmd'%(ts.frame),
    #                      sphpdb='hole2_output/hole%i.sph'%(ts.frame),
    #                     #  sph_process='~/hole2/exe/sph_process',
    #                   	  sos_triangle='~/hole2/exe/sos_triangle')

    # rxn_coords = profiles[0].rxn_coord
    # pore_length = rxn_coords[-1] - rxn_coords[0]
    pore_radius = profiles[0].radius
    pore_list = []
    with open(working_dir+'/radius_PROFILE_%s.out'%(chain),"a") as file_out:
        with open(working_dir+'/hole%i.out'%(frame),"r") as file_one:
            patrn = "sampled"
            for line in file_one:
                if re.search(patrn,line):
                    line = line.split()
                    pore_list = np.append(pore_list,(line[0],line[1]))
                    file_out.write("%s\t%s\t%s\n"%(frame,line[0],line[1]))
                file_out.flush()


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

chainA = "protein and segid PROA"
chainB = "protein and segid PROB"

for ts in u.trajectory[::traj_skip]:
    print(ts.frame)
    pore_cen_A = u.select_atoms("%s"%(chainA)).center_of_geometry()
    pore_cen_B = u.select_atoms("%s"%(chainB)).center_of_geometry()
    pore_radius_A = pore_size(pore_cen_A,ts.frame,"A")
    pore_radius_B = pore_size(pore_cen_B,ts.frame,"B")
