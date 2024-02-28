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

def define_topology(AtomGroup):
    all_min_z=[]
    all_max_z=[]
    for ts in u.trajectory[::traj_skip]:
        all_min_z.append(AtomGroup.positions[:,2].min())
        all_max_z.append(AtomGroup.positions[:,2].max())
    return(np.average(all_min_z),np.average(all_max_z))

def define_new_topology(sel1,sel2,sel3,sel4):
    all_z = [*sel1,*sel2,*sel3,*sel4]
    new_defined_topology = (all_z[0],np.average(all_z[1:3]),np.average(all_z[3:5]),np.average(all_z[5:7]),all_z[7])
    return(np.round(all_z,1))
    
# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--skip', type=int, default=10,
                    help='Skipping rate  Default = 10 frames')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='System Name')

args = parser.parse_args()

top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
systemname = args.system

u = mda.Universe(top_file,traj_file)

ca_binding_chainA = u.select_atoms("protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729 and name CA and segid PROA",updating=True)
inner_vest_chainA = u.select_atoms("protein and resid 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710 and name CA and segid PROA",updating=True)
neck_chainA = u.select_atoms("protein and resid 505 506 507 508 540 54 542 543 590 591 592 593 635 636 711 712 and name CA and segid PROA",updating=True)
outer_vest_chainA = u.select_atoms("protein and resid 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716 and name CA and segid PROA",updating=True)

ca_binding_chainB = u.select_atoms("protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729 and name CA and segid PROB",updating=True)
inner_vest_chainB = u.select_atoms("protein and resid 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710 and name CA and segid PROB",updating=True)
neck_chainB = u.select_atoms("protein and resid 505 506 507 508 540 54 542 543 590 591 592 593 635 636 711 712 and name CA and segid PROB",updating=True)
outer_vest_chainB = u.select_atoms("protein and resid 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716 and name CA and segid PROB",updating=True)


ca_binding_topology_chainA = define_topology(ca_binding_chainA)
inner_vest_topology_chainA = define_topology(inner_vest_chainA)
neck_topology_chainA = define_topology(neck_chainA)
outer_vest_topology_chainA = define_topology(outer_vest_chainA)

ca_binding_topology_chainB = define_topology(ca_binding_chainB)
inner_vest_topology_chainB = define_topology(inner_vest_chainB)
neck_topology_chainB = define_topology(neck_chainB)
outer_vest_topology_chainB = define_topology(outer_vest_chainB)

chaina = define_new_topology(ca_binding_topology_chainA,inner_vest_topology_chainA,neck_topology_chainA,outer_vest_topology_chainA)
chainb = define_new_topology(ca_binding_topology_chainB,inner_vest_topology_chainB,neck_topology_chainB,outer_vest_topology_chainB)

with open("PORE_TOPOLOGY_%s.dat"%(systemname),'w+') as TOPOLOGY:
    # print("Ca_binding","A",ca_binding_topology_chainA)
    TOPOLOGY.write("#NAME CHAIN MIN_Z MAX_Z \n")
    TOPOLOGY.write("%s\t%s\t%s\t%s\n"%("Ca_binding","A",chaina[0],chaina[1]))
    TOPOLOGY.write("%s\t%s\t%s\t%s\n"%("Inner_Vest","A",chaina[1],chaina[2]))
    TOPOLOGY.write("%s\t%s\t%s\t%s\n"%("Neck","A",chaina[2],chaina[3]))
    TOPOLOGY.write("%s\t%s\t%s\t%s\n"%("Outer_Vest","A",chaina[3],chaina[4]))

    TOPOLOGY.write("%s\t%s\t%s\t%s\n"%("Ca_binding","B",chainb[0],chainb[1]))
    TOPOLOGY.write("%s\t%s\t%s\t%s\n"%("Inner_Vest","B",chainb[1],chainb[2]))
    TOPOLOGY.write("%s\t%s\t%s\t%s\n"%("Neck","B",chainb[2],chainb[3]))
    TOPOLOGY.write("%s\t%s\t%s\t%s\n"%("Outer_Vest","B",chainb[3],chainb[4]))
