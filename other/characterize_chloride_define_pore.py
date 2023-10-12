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

sns.set_context("paper")

def find_zone(AtomGroup):
    minZ = AtomGroup.positions[:,2].min()
    maxZ = AtomGroup.positions[:,2].max()
    return minZ,maxZ


# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--index', type=int, nargs='*',
                    help='Indexes of tracked Chloride(s) as a list, separated by a space')

parser.add_argument('--chain', type=str,
                    help='Chain of protein')

parser.add_argument('--pore_traj', type=str,default=None,
                    help='Input trajectory for REGION identifications')

parser.add_argument('--pore_skip', type=int, default=100,
                    help='Skipping rate for REGION identifications. Default = 10 frames')

parser.add_argument('--combined', action='store_true', default=False,
                    help='Skipping rate for REGION identifications. Default = 10 frames')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
index_list = args.index
chain_name = args.chain
pore_traj = args.pore_traj
pore_skip = args.pore_skip
combined = args.combined
# print(top_file,traj_file,index_list,chain_name)


u = mda.Universe(top_file,traj_file)

with open("track_all_cl_at_chain_%s.dat"%(chain_name.upper()),'w+') as of:
    of.write("#frame index z_cl @chain\n")
    for ind in index_list:
        chloride = u.select_atoms("index %s"%(ind),updating=True)
        for ts in u.trajectory[::10]:
            z_cl = chloride.positions[:,2]
            for i in np.arange(len(z_cl)):
                of.write("%s\t%s\t%s\t%s\n"%(ts.frame,ind,z_cl[i],chain_name.upper()))
    print("REMEMBER TO CONCATENATE TWO track_all_cl_at_chain_*.DAT FILES FOR THE NEXT ANALYSIS")

if pore_traj!=None:
    u_pore = mda.Universe(top_file,pore_traj)

    ca_binding = u_pore.select_atoms("protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729 and backbone",updating=True)
    inner_vest = u_pore.select_atoms("protein and resid 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710 725 and backbone",updating=True)
    neck = u_pore.select_atoms("protein and resid 505 506 507 508 540 54 542 543 590 591 592 593 635 636 711 712 and backbone",updating=True)
    outer_vest = u_pore.select_atoms("protein and resid 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716 and backbone",updating=True)
    with open("define_pore.dat",'w+') as of:
        of.write('''#protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729
#protein and 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710 725
#protein and 505 506 507 508 540 54 542 543 590 591 592 593 635 636 711 712
#protein and 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716\n''')
        of.write("# Frame ca_binding_min,ca_binding_max,inner_vest_min,inner_vest_max,neck_min,neck_max,outer_vest_min,outer_vest_max\n")
        for ts in u_pore.trajectory[::pore_skip]:
            ca_binding_min,ca_binding_max = find_zone(ca_binding)
            inner_vest_min,inner_vest_max = find_zone(inner_vest)
            neck_min,neck_max             = find_zone(neck)
            outer_vest_min,outer_vest_max = find_zone(outer_vest)
            # print(ca_binding_min,ca_binding_max,inner_vest_min,inner_vest_max,neck_min,neck_max,outer_vest_min,outer_vest_max)
            of.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(ts.frame,ca_binding_min,ca_binding_max,inner_vest_min,inner_vest_max,neck_min,neck_max,outer_vest_min,outer_vest_max))
else:
    print("NO PORE DEFINITION IS CALCULATED!")