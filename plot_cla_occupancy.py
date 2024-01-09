import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc

sns.set(style="ticks", context="notebook")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to Rg profile')

parser.add_argument('--cutoff', type=int, default='50',
                    help='Cutoff distance from the pore region to search for Cl-. Default 50 (A)')
args = parser.parse_args()

ifile =  args.input
cutoff = args.cutoff

data = pd.read_csv(ifile,comment='#',
delim_whitespace=True,
names=['FRAME','INDEX','CHAIN','ZPOS','POSITION_PORE','DCOM_PORE'])
data['ns']=data.FRAME*0.2
data['COUNT']=1
data_count = data.groupby(['CHAIN','POSITION_PORE']).sum().reset_index()
data_count2 = data_count[['CHAIN','POSITION_PORE','COUNT']].query("POSITION_PORE =='Inner_Vestibule' or POSITION_PORE =='Ca_Binding'")
data_count3=data_count2.groupby(['CHAIN']).sum().reset_index()
# print(data_count3)
if len(data_count3.CHAIN.unique())==1:
    if data_count3.CHAIN.unique()=='A':
        # print("Mop")
        chain_dummy={'CHAIN':['B'],'COUNT':0}
        data_count3=pd.concat([data_count3,pd.DataFrame(chain_dummy)],ignore_index=True)
    else:
        chain_dummy={'CHAIN':['A'],'COUNT':0}
        data_count3=pd.concat([data_count3,pd.DataFrame(chain_dummy)],ignore_index=True)
# print(data_count3)
data_count3['Occupancy']=data_count3['COUNT']/10000*100
g = sns.barplot(data=data_count3.sort_values('CHAIN'),
x='CHAIN',
y='Occupancy',
# hue='CHAIN',
)
weight='bold'
fontsize=30
plt.xlabel("Z Coordinate")
plt.xticks(size=fontsize,weight=weight)
plt.yticks(size=fontsize,weight=weight)
plt.ylabel("Occupancy %",size=fontsize,weight=weight)
plt.ylim([0,70])

# # # ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,10)

plt.tight_layout()
plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
# plt.show()
