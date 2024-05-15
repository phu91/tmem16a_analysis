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
data['COUNT']=1
data_count_a = data[['FRAME','CHAIN','POSITION_PORE','COUNT']].query("POSITION_PORE =='Inner_Vestibule' or POSITION_PORE =='Ca_Binding'")
# print(data_count_a)
n_frame = len(data['FRAME'].unique())
OCCUPANCY=round(len(data_count_a)/n_frame*100,2)
# print(OCCUPANCY)
print("\nFile: %s"%(ifile))
print("Number of FRAME: %s"%(n_frame))
print("Percent of Finding Cl-: %s\n"%(OCCUPANCY))

### UNCOMMENT THIS SECTION TO GET THE OCCUPANCY FOR EACH CHAIN
# data_count_b = data_count_a.groupby(['CHAIN']).sum()
# # print(data_count_b)
# n_frame = len(data['FRAME'].unique())
# data_count_b['OCCUPANCY']=data_count_b['COUNT']/n_frame*100
# # print(data_count_b)
# data_count_b = data_count_b.reset_index()
# data_count_c = data_count_b[['CHAIN','OCCUPANCY']]
# print(data_count_c)


# g = sns.barplot(data=data_count_c.sort_values('CHAIN'),
# x='CHAIN',
# y='OCCUPANCY',
# # hue='CHAIN',
# )
# weight='bold'
# fontsize=30
# plt.xlabel("Z Coordinate")
# plt.xticks(size=fontsize,weight=weight)
# plt.yticks(size=fontsize,weight=weight)
# plt.ylabel("Occupancy %",size=fontsize,weight=weight)
# plt.ylim([0,70])

# # # # ### MISCELLANEOUS ###
# plt.suptitle("%s"%(ifile[:-4]),va='top')
# plt.rcParams['ps.useafm'] = True
# # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,10)

# plt.tight_layout()
# plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
# plt.show()
