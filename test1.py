import MDAnalysis as mda
import pandas as pd
import numpy as np
import argparse
import warnings
# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')


# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--input', type=str,
                    help='A required topology (PSF, PDB, etc.) file')
parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Add a system name to output file')
args = parser.parse_args()


ifile =  args.input
systemname = args.system


data2 = pd.read_csv(ifile,
names=['frame','resid','position','binding'],sep='\t')

# print(data2)
residue_count=len(data2)
residue_id=data2.resid
residue_position=data2.position
# print(residue_id)
# print(residue_position)
print(residue_count)
# FILTERED DATA
filtered_df = pd.DataFrame()
for i in range(residue_count):
    filtered = data2.query("resid==%s and position=='%s'"%(residue_id[i],residue_position[i]))
    # print(filtered)
    filtered.loc[filtered['position']=='567B_887A','position']="TM10_A"
    filtered.loc[filtered['position']=='567A_887B','position']="TM10_B"
    filtered_df = pd.concat([filtered,filtered_df],ignore_index=True)
# print(filtered_df)
filtered_df = filtered_df.sort_values(by=['frame','position'])
# print(filtered_df)
u=filtered_df.groupby(['frame','position']).sum().reset_index()
u = u[['frame','position','binding']]
# print(u)
u.loc[u['binding']>=1,'binding']=1  ## If any PIP2 is at this position, BINDING=1 as sum(binding)>1
# print(u)
u.to_csv("PIP2_TM345_10_PROFILE_%s_SIMPLE.dat"%(systemname),index=False,sep='\t',header=False,chunksize=5000)

# with open("PIP2_TM345_10_PROFILE_%s_SIMPLE.dat"%(systemname),'r+') as ofile:
#     lines = ofile.readlines()     # lines is list of line, each element '...\n'
#     # lines.insert(0, one_line)  # you can use any index if you know the line index
#     ofile.seek(0)                 # file pointer locates at the beginning to write the whole file again
#     ofile.write("## PIP2 FOUND BETWEEN K567(LOOP45) to K887(THIRD_CA). CUTOFF: %s A\n"%(cutoff))
#     ofile.write("## FOUND PIP2 LIST:\n")
#     ofile.write("## resid position\n")
#     # for ind,row in df_pip2.iterrows():
#         # ofile.write("#\t%s\t%s\n"%(row.resid,row.position))
#     ofile.write("# FRAME RESID POSITION BINDING\n")
#     ofile.writelines(lines)       # write whole lists again to the same file