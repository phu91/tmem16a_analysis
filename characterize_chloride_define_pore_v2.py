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

parser.add_argument('--skip', type=int, default=10,
                    help='Skipping rate  Default = 10 frames')

parser.add_argument('--pore',action='store_true',default=False,
                    help='Calculate PORE TOPOLOGY')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='System Name')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
pore_traj = args.pore
systemname = args.system

u = mda.Universe(top_file,traj_file)

all = u.select_atoms("all")

chloride_list = []
with mda.Writer('cla_pore_a.dcd', all.n_atoms) as w:
    for ts in u.trajectory[::traj_skip]:
        print(ts.frame)
        chloride_a = u.select_atoms("(around 5 (resid 550 551 584 641 645 646 700 701 and protein and segid PROA)) and resname CLA",updating=True)
        if len(chloride_a)!=0:
            for ion in chloride_a.atoms.indices:
                chloride_list.append((ion,"A"))
                w.write(all)

print()

with mda.Writer('cla_pore_b.dcd', all.n_atoms) as w:
    for ts in u.trajectory[::traj_skip]:
        print(ts.frame)
        chloride_b = u.select_atoms("(around 5 (resid 550 551 584 641 645 646 700 701 and protein and segid PROB)) and resname CLA",updating=True)
        if len(chloride_b)!=0:
            for ion in chloride_b.atoms.indices:
                chloride_list.append((ion,"B"))
                w.write(all)

# print(len(chloride_list))
df_cla = pd.DataFrame(chloride_list,columns=['atom_index', 'chain']).drop_duplicates()

print(df_cla)
with open("cla_track_all_%s.dat"%(systemname),'w+') as of:
    of.write("#frame index at_chain z_pos distance_com\n")
    for row in df_cla.itertuples():
        # print(row[2])
        com_pore_a = u.select_atoms("resid 550 551 584 641 645 646 700 701 and protein and segid PROA", updating=True).center_of_mass()
        com_pore_b = u.select_atoms("resid 550 551 584 641 645 646 700 701 and protein and segid PROB", updating=True).center_of_mass()

        selected_chloride = u.select_atoms("index %s"%(row[1]),updating=True)
        for ts in u.trajectory[::traj_skip]:
            # print(ts.frame,*selected_chloride.atoms.indices,row[2],selected_chloride.positions[0,2])
            if row[2] == "A":
                of.write("%s\t%s\t%s\t%s\t%s\n"%(ts.frame,*selected_chloride.atoms.indices,"A",selected_chloride.positions[0,2],abs(selected_chloride.positions[0,2]-com_pore_a[2])))
                of.flush()
            else:
                of.write("%s\t%s\t%s\t%s\t%s\n"%(ts.frame,*selected_chloride.atoms.indices,"B",selected_chloride.positions[0,2],abs(selected_chloride.positions[0,2]-com_pore_b[2])))
                of.flush()
if pore_traj is True:
    print("\nPORE TOPOLOGY IS CALCULATING!\n")
    ca_binding = u.select_atoms("protein and resid 643 648 649 651 652 653 654 700 702 703 704 705 727 729 and backbone",updating=True)
    inner_vest = u.select_atoms("protein and resid 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710 and backbone",updating=True)
    neck = u.select_atoms("protein and resid 505 506 507 508 540 541 542 543 590 591 592 593 635 636 711 712 and backbone",updating=True)
    outer_vest = u.select_atoms("protein and resid 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716 and backbone",updating=True)
    with open("define_pore_%s.dat"%(systemname),'w+') as of:
        of.write('''#protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729
#protein and 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710
#protein and 505 506 507 508 540 54 542 543 590 591 592 593 635 636 711 712
#protein and 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716\n''')
        of.write("# Frame ca_binding_min,ca_binding_max,inner_vest_min,inner_vest_max,neck_min,neck_max,outer_vest_min,outer_vest_max\n")
        for ts in u.trajectory[::traj_skip]:
            ca_binding_min,ca_binding_max = find_zone(ca_binding)
            inner_vest_min,inner_vest_max = find_zone(inner_vest)
            neck_min,neck_max             = find_zone(neck)
            outer_vest_min,outer_vest_max = find_zone(outer_vest)
            # print(ca_binding_min,ca_binding_max,inner_vest_min,inner_vest_max,neck_min,neck_max,outer_vest_min,outer_vest_max)
            of.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(ts.frame,ca_binding_min,ca_binding_max,inner_vest_min,inner_vest_max,neck_min,neck_max,outer_vest_min,outer_vest_max))
            of.flush()
else:
    print("\nPORE TOPOLOGY IS NOT CALCULATED!\n")
