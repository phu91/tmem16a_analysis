# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader
import math, sys
import argparse, math
from tqdm import tqdm

# import multiprocessing
# from multiprocessing import Pool
# from functools import partial

def find_zone(AtomGroup):
    minZ = AtomGroup.positions[:,2].min()
    maxZ = AtomGroup.positions[:,2].max()
    return minZ,maxZ

def position_name(ca_binding_current,inner_vest_current,neck_current,outer_vest_current,chloride_pos_current):
    position='Out_of_Zone'
    if ca_binding_current[0]<=chloride_pos_current[2]<=ca_binding_current[1]:
        # print(ts.frame,chloride_pos_current[2],ca_binding)
        position='Ca_Binding'
    elif inner_vest_current[0]<=chloride_pos_current[2]<=inner_vest_current[1]:
        # print(ts.frame,chloride_pos_current[2],inner_vest)
        position='Inner_Vestibule'
    elif neck_current[0]<=chloride_pos_current[2]<=neck_current[1]:
        # print(ts.frame,chloride_pos_current[2],neck)
        position='Neck'
    elif outer_vest_current[0]<=chloride_pos_current[2]<=outer_vest_current[1]:
        # print(ts.frame,chloride_pos_current[2],outer_vest)
        position='Outer_Vestibule'
    return position

def find_chlorides():
    # print("\n===> CUT-OFF distance of CHLORIDES and the protein is %s Angstroms\n"%(cutoff))
    chloride_list_found = []
    with mda.Writer('CHLORIDES_at_pore_chain_A_%s.dcd'%(systemname), all.n_atoms) as w:
        for ts in tqdm(u.trajectory[::traj_skip],desc='Finding Cl- at chain A'):
            # print(ts.frame)
            chloride_a = u.select_atoms("(around %s (resid 550 551 584 641 645 646 700 701 and segid PROA)) and resname CLA"%(cutoff),updating=True)
            if len(chloride_a)!=0:
                # print(chloride_a)
                for ion in chloride_a.atoms.indices:
                    chloride_list_found.append((ion,"A"))
                    w.write(all)

    with mda.Writer('CHLORIDES_at_pore_chain_B_%s.dcd'%(systemname), all.n_atoms) as w:
        for ts in tqdm(u.trajectory[::traj_skip],desc='Finding Cl- at chain B'):
            # print(ts.frame)
            chloride_b = u.select_atoms("(around %s (resid 550 551 584 641 645 646 700 701 and segid PROB)) and resname CLA"%(cutoff),updating=True)
            if len(chloride_b)!=0:
                for ion in chloride_b.atoms.indices:
                    chloride_list_found.append((ion,"B"))
                    w.write(all)
    print("\n===> Trajectories of FOUND Cl- are saved in 'Chlorides_at_pore_chain_A_%s.dcd' and 'Chlorides_at_pore_chain_B_%s.dcd'\n"%(systemname,systemname))
    return chloride_list_found

# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate  Default = 1 frames')

parser.add_argument('--cutoff', type=int, default='5',
                    help='Cut-off distance to find Cl-. Default is 5A')
                    
parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='System Name')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
cutoff = args.cutoff
systemname = args.system

u = mda.Universe(top_file,traj_file)

all = u.select_atoms("all")

### FIND CHLORIDES 

chloride_list = find_chlorides()
# print(chloride_list)
df_chlorides = pd.DataFrame(chloride_list,columns=['atom_index','chain']).drop_duplicates()
print("\n===> Found %s chlorides with a cut-off of %s Angstrom"%(len(df_chlorides),cutoff))
print(df_chlorides.reset_index())
### TRACK EACH CHLRODIE

with open("CHLORIDES_PROFILE_%s.dat"%(systemname),'w+') as of:
    of.write("#frame index at_chain z_coord at_position distance_com_pore\n")
    for row in df_chlorides.itertuples():
        # print(row[2])
        com_pore_a = u.select_atoms("resid 550 551 584 641 645 646 700 701 and segid PROA", updating=True).center_of_mass()
        com_pore_b = u.select_atoms("resid 550 551 584 641 645 646 700 701 and segid PROB", updating=True).center_of_mass()
        selected_chloride = u.select_atoms("index %s"%(row[1]),updating=True)
        # print(selected_chloride)
        for ts in u.trajectory[::traj_skip]:
            ca_binding = u.select_atoms("protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729 and backbone",updating=True)
            inner_vest = u.select_atoms("protein and resid 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710 and backbone",updating=True)
            neck = u.select_atoms("protein and resid 505 506 507 508 540 54 542 543 590 591 592 593 635 636 711 712 and backbone",updating=True)
            outer_vest = u.select_atoms("protein and resid 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716 and backbone",updating=True)
            ca_binding = find_zone(ca_binding)
            inner_vest = find_zone(inner_vest)
            neck       = find_zone(neck)
            outer_vest = find_zone(outer_vest)
            # print(ts.frame,*selected_chloride.atoms.indices,row[2],selected_chloride.positions)

            chloride_coords=selected_chloride.positions[0]
            current_position = position_name(ca_binding,inner_vest,neck,outer_vest,chloride_coords)

            if row[2] == "A":
                    of.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(ts.frame,*selected_chloride.atoms.indices,"A",chloride_coords[2],current_position,math.dist(chloride_coords,com_pore_a)))
                    of.flush()
            else:
                of.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(ts.frame,*selected_chloride.atoms.indices,"B",chloride_coords[2],current_position,math.dist(chloride_coords,com_pore_b)))
                of.flush()
    print("\n===> Output 'CHLORIDES_PROFILE_%s.dat' is saved\n"%(systemname))
