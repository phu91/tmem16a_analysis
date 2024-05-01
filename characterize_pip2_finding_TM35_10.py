import MDAnalysis as mda
import pandas as pd
import numpy as np
import math, sys
import argparse
from MDAnalysis.coordinates.memory import MemoryReader
from tqdm import tqdm
from MDAnalysis.analysis import dihedrals
import warnings
from MDAnalysis.analysis import distances

# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

def distance_by_chain(frameNow,sel1,sel2,chainName):
    oneDistance = distances.distance_array(sel1.positions, # reference
                                    sel2.positions, # configuration
                                    box=u.dimensions)
    dict_oneDistance = {'frame':frameNow,
                        'resid':sel1.resids,
                        'chain':chainName,
                        'distance':np.round(oneDistance[0],2)
                        }

    df_oneDistance = pd.DataFrame(dict_oneDistance)
    return df_oneDistance



# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, PDB, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--begin', type=int, default=1,
                    help='Starting Frame. Default = 1 FRAME 1')

parser.add_argument('--end', type=int, default=-1,
                    help='Ending Frame. Default = -1 ALL FRAME')

parser.add_argument('--cutoff', type=float, default=4,
                    help='Cutoff. Default = 4 A')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate for Radial Distribution Function calculations. Default = 1 frames')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Add a system name to output file')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
systemname = args.system
cutoff = args.cutoff

print("LOAD DATA ...\n")
u = mda.Universe(top_file,traj_file)

K567_A = u.select_atoms("segid PROA and resid 567 and name NZ",updating=True)
K887_A = u.select_atoms("segid PROA and resid 887 and name NZ",updating=True)
K567_B = u.select_atoms("segid PROB and resid 567 and name NZ",updating=True)
K887_B = u.select_atoms("segid PROB and resid 887 and name NZ",updating=True)

# important_pip2_A = u.select_atoms("(around 5 ((segid PROA and resid 567 and name NZ) and (segid PROB and resid 887 and name NZ)) and (resname PLPI24)",updating=True)
# important_pip2_B = u.select_atoms("(around 5 ((segid PROB and resid 567 and name NZ) and (segid PROA and resid 887 and name NZ)) and (resname PLPI24)",updating=True)
important_pip2_A = u.select_atoms("(around %s ((segid PROA and resid 567 and name NZ) or (segid PROB and resid 887 and name NZ))) and (resname PLPI24 and name P4)"%(cutoff),updating=True)
important_pip2_B = u.select_atoms("(around %s ((segid PROB and resid 567 and name NZ) or (segid PROA and resid 887 and name NZ))) and (resname PLPI24 and name P4)"%(cutoff),updating=True)

df = pd.DataFrame()

pip2_list_found =[]
for ts in tqdm(u.trajectory[::traj_skip*5],desc='Finding PIP2'):
        if len(important_pip2_A)!=0:
            for pip2 in important_pip2_A.resids:
                pip2_list_found.append((pip2,"567A_887B"))
        if len(important_pip2_B)!=0:
            for pip2 in important_pip2_B.resids:
                pip2_list_found.append((pip2,"567B_887A"))
# print(pip2_list_found)

df_pip2 = pd.DataFrame(pip2_list_found,columns=['resid','position']).drop_duplicates()
if len(df_pip2)!=0:
    print()
    print(df_pip2)
else:
    "No PIP2 found.\n"
print()

with open("PIP2_TM345_10_PROFILE_%s.dat"%(systemname),'w+') as ofile:
    if len(df_pip2)!=0:
        for ts in tqdm(u.trajectory[::traj_skip],desc='Tracking PIP2'):
            for res in df_pip2.resid:
                pip2_A = u.select_atoms("(around 10 ((segid PROA and resid 567 and name NZ) or (segid PROB and resid 887 and name NZ))) and (resname PLPI24 and name P4 and resid %s)"%(res),updating=True)
                pip2_B = u.select_atoms("(around 10 ((segid PROB and resid 567 and name NZ) or (segid PROA and resid 887 and name NZ))) and (resname PLPI24 and name P4 and resid %s)"%(res),updating=True)
                # print(pip2_A.resids)
                # print(pip2_B.resids)
                if len(pip2_A.resids)!=0:
                    # print(ts.frame,res,"567A_887B","1")
                    ofile.write("%s\t%s\t%s\t%s\n"%(ts.frame,res,"567A_887B","1"))
                else:
                    # print(ts.frame,res,"567A_887B","0")
                    ofile.write("%s\t%s\t%s\t%s\n"%(ts.frame,res,"567A_887B","0"))
                if len(pip2_B.resids)!=0:
                    # print(ts.frame,res,"567B_887A","1")
                    ofile.write("%s\t%s\t%s\t%s\n"%(ts.frame,res,"567B_887A","1"))

                else:
                    # print(ts.frame,res,"567B_887A","0")
                    ofile.write("%s\t%s\t%s\t%s\n"%(ts.frame,res,"567B_887A","0"))
    else:
        print("No PIP2 found.\n")


data2 = pd.read_csv("PIP2_TM345_10_PROFILE_%s.dat"%(systemname),
names=['frame','resid','position','binding'],sep='\t')

residue_count=len(df_pip2)
residue_id=df_pip2.resid
residue_position=df_pip2.position
# print(data2)
## FILTERED DATA
filtered_df = pd.DataFrame()
for i in range(residue_count):
    filtered = data2.query("resid==%s and position=='%s'"%(residue_id[i],residue_position[i]))
    filtered.loc[filtered['position']=='567B_887A','position']="TM10_A"
    filtered.loc[filtered['position']=='567A_887B','position']="TM10_B"
    filtered_df = pd.concat([filtered,filtered_df],ignore_index=True)
# print(filtered_df)
filtered_df = filtered_df.sort_values(by=['frame','position'])
print(filtered_df)
filtered_df.to_csv("PIP2_TM345_10_PROFILE_%s.dat"%(systemname),index=False,sep='\t',header=False)

with open("PIP2_TM345_10_PROFILE_%s.dat"%(systemname),'r+') as ofile:
    lines = ofile.readlines()     # lines is list of line, each element '...\n'
    # lines.insert(0, one_line)  # you can use any index if you know the line index
    ofile.seek(0)                 # file pointer locates at the beginning to write the whole file again
    ofile.write("## PIP2 FOUND BETWEEN K567(LOOP45) to K887(THIRD_CA). CUTOFF: %s A\n"%(cutoff))
    ofile.write("## FOUND PIP2 LIST:\n")
    ofile.write("## resid position\n")
    for ind,row in df_pip2.iterrows():
        ofile.write("#\t%s\t%s\n"%(row.resid,row.position))
    ofile.write("# FRAME RESID POSITION BINDING\n")
    ofile.writelines(lines)       # write whole lists again to the same file

# with open("PIP2_TM345_10_PROFILE_%s.dat"%(systemname),'w+') as ofile:
    # ofile.write("## PIP2 FOUND BETWEEN K567(LOOP45) to K887(THIRD_CA). CUTOFF: %s A\n"%(cutoff))
    # ofile.write("## FOUND PIP2 LIST. Use this list to filter out the unrelevant rows.\n")
    # ofile.write("# resid position\n")
    # for ind,row in df_pip2.iterrows():
    #     ofile.write("#\t%s\t%s\n"%(row.resid,row.position))
    # ofile.write("# FRAME RESID POSITION BINDING\n")