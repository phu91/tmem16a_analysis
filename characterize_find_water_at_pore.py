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


# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate  Default = 1 frames')

parser.add_argument('--cutoff', type=float, default='3.5',
                    help='Cut-off distance to find water. Default is 3.5A')
                    
parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='System Name')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
cutoff = args.cutoff
systemname = args.system

u = mda.Universe(top_file,traj_file)

water_chainA_ca_binding = u.select_atoms("name OH2 and (around %s (protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729 and backbone and segid PROA))"%(cutoff),updating=True)
water_chainA_inner_vest = u.select_atoms("name OH2 and (around %s (protein and resid 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710 and backbone and segid PROA))"%(cutoff),updating=True)
water_chainA_neck       = u.select_atoms("name OH2 and (around %s (protein and resid 505 506 507 508 540 542 543 590 591 592 593 635 636 711 712 and backbone and segid PROA))"%(cutoff),updating=True)
water_chainA_outer_vest = u.select_atoms("name OH2 and (around %s (protein and resid 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716 and backbone and segid PROA))"%(cutoff),updating=True)
water_chainB_ca_binding = u.select_atoms("name OH2 and (around %s (protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729 and backbone and segid PROB))"%(cutoff),updating=True)
water_chainB_inner_vest = u.select_atoms("name OH2 and (around %s (protein and resid 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710 and backbone and segid PROB))"%(cutoff),updating=True)
water_chainB_neck       = u.select_atoms("name OH2 and (around %s (protein and resid 505 506 507 508 540 542 543 590 591 592 593 635 636 711 712 and backbone and segid PROB))"%(cutoff),updating=True)
water_chainB_outer_vest = u.select_atoms("name OH2 and (around %s (protein and resid 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716 and backbone and segid PROB))"%(cutoff),updating=True)

# print(water_chainA_ca_binding.resids)
# print(water_chainA_inner_vest.resids)
df = pd.DataFrame()

for ts in tqdm(u.trajectory[::traj_skip]):
#### CHAIN A
    data_1 = {'frame':ts.frame,
            'resid':water_chainA_ca_binding.resids,
            'position':'CA_BINDING',
            'chain':'A'
            }
    data_2 = {'frame':ts.frame,
            'resid':water_chainA_inner_vest.resids,
            'position':'INNER_VEST',
            'chain':'A'
            }
    data_3 = {'frame':ts.frame,
            'resid':water_chainA_neck.resids,
            'position':'NECK',
            'chain':'A'
            }
    data_4 = {'frame':ts.frame,
            'resid':water_chainA_outer_vest.resids,
            'position':'OUTER_VEST',
            'chain':'A'
            }
#### CHAIN B
    data_5 = {'frame':ts.frame,
            'resid':water_chainB_ca_binding.resids,
            'position':'CA_BINDING',
            'chain':'B'
            }
    data_6 = {'frame':ts.frame,
            'resid':water_chainB_inner_vest.resids,
            'position':'INNER_VEST',
            'chain':'B'
            }
    data_7 = {'frame':ts.frame,
            'resid':water_chainB_neck.resids,
            'position':'NECK',
            'chain':'B'
            }
    data_8 = {'frame':ts.frame,
            'resid':water_chainB_outer_vest.resids,
            'position':'OUTER_VEST',
            'chain':'B'
            }

    data_1_A_df = pd.DataFrame(data_1)
    data_2_A_df = pd.DataFrame(data_2)
    data_3_A_df = pd.DataFrame(data_3)
    data_4_A_df = pd.DataFrame(data_4)
    data_5_B_df = pd.DataFrame(data_5)
    data_6_B_df = pd.DataFrame(data_6)
    data_7_B_df = pd.DataFrame(data_7)
    data_8_B_df = pd.DataFrame(data_8)
    # print(data_1_df)
    df = pd.concat([df,data_1_A_df,data_2_A_df,data_3_A_df,data_4_A_df,data_5_B_df,data_6_B_df,data_7_B_df,data_8_B_df])

df = df.drop_duplicates(['resid','position'])
df = df.reset_index().drop(columns=['index','resid'])

df['count']=1
df = df.groupby(['frame','position','chain']).sum()
df = df.reset_index()

df.to_csv("HYDRATION_PROFILE_%s.dat"%(systemname),sep='\t',index=False)
print("\n==> DATA saved to HYDRATION_PROFILE_%s.dat" %(systemname))