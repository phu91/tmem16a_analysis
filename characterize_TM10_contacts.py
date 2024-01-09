# Updated on 11/29/2022
import MDAnalysis as mda
import numpy as np
import math, sys
import argparse
from MDAnalysis.analysis import align, distances
from MDAnalysis.coordinates.memory import MemoryReader
from tqdm import tqdm
import pandas as pd

def where_is_this_residue(RESID):
    RESID=int(RESID)
    ## Define Helix
    N_term = np.arange(0,327)
    helix_01_range = np.arange(327,361)
    helix_02_range = np.arange(399,440)
    helix_03_range = np.arange(478,519)
    helix_04_range = np.arange(534,563)
    helix_05_range = np.arange(572,599)
    helix_06_range = np.arange(630,667)
    helix_07_range = np.arange(692,714)
    helix_08_range = np.arange(718,741)
    helix_09_range = np.arange(753,781)
    helix_10_range_cterm = np.arange(887,928)

    # print(helix_01_range)

    ## Define Loop
    loop_1_2_range  = np.arange(361,399)
    loop_2_3_range  = np.arange(440,478)
    loop_3_4_range  = np.arange(519,534)
    loop_4_5_range  = np.arange(563,572)
    loop_5_6_range  = np.arange(599,630)
    loop_6_7_range  = np.arange(667,692)
    loop_7_8_range  = np.arange(714,718)
    loop_8_9_range  = np.arange(741,753)
    loop_9_10_range_tm10_top = np.arange(781,887)
    PART_NAME='UNK'
    if RESID in helix_01_range:
        PART_NAME='TM1'
    elif RESID in helix_02_range:
        PART_NAME='TM2'
    elif RESID in helix_03_range:
        PART_NAME='TM3'
    elif RESID in helix_04_range:
        PART_NAME='TM4'
    elif RESID in helix_05_range:
        PART_NAME='TM5'
    elif RESID in helix_06_range:
        PART_NAME='TM6'
    elif RESID in helix_07_range:
        PART_NAME='TM7'
    elif RESID in helix_08_range:
        PART_NAME='TM8'
    elif RESID in helix_09_range:
        PART_NAME='TM9'
    elif RESID in helix_10_range_cterm:
        PART_NAME='TM10_CTERM'
    elif RESID in loop_1_2_range:
        PART_NAME='L12'
    elif RESID in loop_2_3_range:
        PART_NAME='L23'
    elif RESID in loop_3_4_range:
        PART_NAME='L34'
    elif RESID in loop_4_5_range:
        PART_NAME='L45'
    elif RESID in loop_5_6_range:
        PART_NAME='L56'
    elif RESID in loop_6_7_range:
        PART_NAME='L67'
    elif RESID in loop_7_8_range:
        PART_NAME='L78'
    elif RESID in loop_8_9_range:
        PART_NAME='L89'
    elif RESID in loop_9_10_range_tm10_top:
        PART_NAME='L910_TM10_TOP'
    elif RESID in N_term:
        PART_NAME='Nterm'
    return PART_NAME

def distance_per_residue(currentFrame,selection_tm10,chain_tm10,selection_helix):

    selection_tm10_positions = selection_tm10.positions
    selection_tm10_atoms = selection_tm10.atoms

    selection_helix_positions = selection_helix.positions
    selection_helix_atoms = selection_helix.atoms

    selection_tm10_label=[]
    for resn,resi in zip(selection_tm10_atoms.resnames,selection_tm10_atoms.resids):
        selection_tm10_label.append((resn+str(resi)))
    # print(len(selection_tm10_label))
    selection_helix_label=[]
    for resn,resi in zip(selection_helix_atoms.resnames,selection_helix_atoms.resids):
        selection_helix_label.append((resn+str(resi)))
    # print(selection_helix_label)

    ## Calculate distance between tm10 and every atoms on the other chain
    dist_arr = distances.distance_array(selection_helix_positions,
                                        selection_tm10_positions,
                                    box=u.dimensions)
    # print(len(selection_tm10_atoms),len(selection_helix_atoms),np.shape(dist_arr))

    ## Add labels for TM10 and selected helix
    distance_df = pd.DataFrame(dist_arr,columns=selection_tm10_label)
    distance_df['OTHER']=selection_helix_label
    distance_df_melt = distance_df.melt(['OTHER'])
    distance_df_melt=distance_df_melt.rename(columns={'variable': "TM10",'value': "DISTANCE"})
    # distance_df_melt_group = distance_df_melt.groupby(['TM10'],sort=False).min()
    # print(distance_df_melt_group)
    ## Select only the cutoff distance, find the smallest distance among the atoms of each residue on OTHER CHAIN
    cut_off_df = distance_df_melt.loc[distance_df_melt.DISTANCE<=cutoff]
    distance_df_melt_cutoff = cut_off_df.groupby(['OTHER'],sort=False).min().round(2)
    distance_df_melt_cutoff['FRAME']=currentFrame
    distance_df_melt_cutoff=distance_df_melt_cutoff.reset_index()
    OTHER_resid = distance_df_melt_cutoff['OTHER'].str[3:]
    # print(OTHER_resid)
    OTHER_PART = OTHER_resid.map(lambda x: where_is_this_residue(x))
    distance_df_melt_cutoff['OTHER_PART']=OTHER_PART
    distance_df_melt_cutoff['TM10_CHAIN']=chain_tm10
    distance_df_melt_cutoff=distance_df_melt_cutoff[['FRAME','TM10_CHAIN','TM10','OTHER','OTHER_PART','DISTANCE']]
    return distance_df_melt_cutoff

# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate for Radial Distribution Function calculations. Default = 1 frames')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Add a system name to output file')

parser.add_argument('--cutoff', type=int,default='5',
                    help='Cut-off distance. Default 5A')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
cutoff = args.cutoff
systemname = args.system

u = mda.Universe(top_file,traj_file,in_memory=True)

# Define NTD of TM10

helix_10a = u.select_atoms('segid PROA and resid 887 to 927 and not name H*',updating=True)
helix_10b = u.select_atoms('segid PROB and resid 887 to 927 and not name H*',updating=True)

chaina = u.select_atoms("segid PROA and not name H*",updating=True)
chainb = u.select_atoms("segid PROB and not name H*",updating=True)

with open ("DISTANCE_TM10_PROFILE_%s"%(systemname),'w+') as ofile:
    ofile.write("FRAME TM10_CHAIN TM10 OTHER OTHER_PART DISTANCE\n")
    collect=[]
    for ts in u.trajectory[::traj_skip]:
        distance_profile_a = distance_per_residue(ts.frame,helix_10a,"A",chainb)
        distance_profile_b = distance_per_residue(ts.frame,helix_10b,"B",chaina)
        collect.append((distance_profile_a.values.tolist()))
        collect.append((distance_profile_b.values.tolist()))

    for i in collect:
        for j in i:
            ofile.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(j[0],j[1],j[2],j[3],j[4],j[5]))
            ofile.flush()