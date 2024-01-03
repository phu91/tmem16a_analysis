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
data_count = data.groupby(['INDEX','CHAIN','POSITION_PORE']).sum().reset_index()
data_count2 = data_count[['INDEX','CHAIN','POSITION_PORE','COUNT']].query("POSITION_PORE !='Out_of_Zone'")
# print(data_count2)
# g = sns.catplot(data=data_count,
# x='POSITION_PORE',
# y='COUNT',
# hue='CHAIN',
# )

barWidth = 0.25
fig = plt.subplots(figsize =(12, 8)) 

species = ("Ca2+ Binding", "Inner Vestibule", "Neck")
# set height of bar 
A = data_count.loc[data_count.CHAIN=='A'].COUNT.tolist()
B = data_count.loc[data_count.CHAIN=='B'].COUNT.tolist()
 
# Set position of bar on X axis 
br1 = np.arange(len(species)) 
br1a = [x + barWidth for x in br1] 
br2 = np.arange(len(species)) 
br2a = [x + barWidth for x in br2] 
# Make the plot
plt.bar(br1a, A, color ='r', width = barWidth, 
        edgecolor ='grey', label ='IT') 
plt.bar(br2a, B, color ='g', width = barWidth, 
        edgecolor ='grey', label ='ECE') 

# Adding Xticks 
# plt.xlabel('Branch', fontweight ='bold', fontsize = 15) 
# plt.ylabel('Students passed', fontweight ='bold', fontsize = 15) 
# plt.xticks([r + barWidth for r in range(len(IT))], 
#         ['2015', '2016', '2017', '2018', '2019'])

### DEFINE PORE TOPOLOGY
# PORE = data[['POSITION_PORE','ZPOS']]
# Inner_Vestibule=PORE.loc[PORE.POSITION_PORE=="Inner_Vestibule"]
# Ca_Binding=PORE.loc[PORE.POSITION_PORE=="Ca_Binding"]
# Neck=PORE.loc[PORE.POSITION_PORE=="Neck"]
# Outer_Vestibule=PORE.loc[PORE.POSITION_PORE=="Outer_Vestibule"]
# Out_of_Zone=PORE.loc[PORE.POSITION_PORE=="Out_of_Zone"]

# # # ### MISCELLANEOUS ###
# # plt.ylim([-30,30])
# plt.suptitle("%s"%(ifile[:-4]),va='top')
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,5.5)
# # plt.locator_params(axis='y', nbins=6)
# plt.legend(loc='upper left',bbox_to_anchor=(1.01, 1),borderaxespad=0)
# plt.tight_layout()
# plt.savefig("%s.png"%(ifile[:-4]))
plt.show()
