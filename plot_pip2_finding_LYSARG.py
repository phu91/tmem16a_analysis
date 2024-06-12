import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from statannot import add_stat_annotation

sns.set_context("notebook")
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT')

args = parser.parse_args()

ifile =  args.input

data = pd.read_csv(ifile,
                 names=['FRAME','CHAIN','RES','ID','COUNT'],
                 comment='#',
                 delim_whitespace=True)
print(data)

# data['TIME']=np.round(data['FRAME']*0.2,1)  ## CHANGE FRAME TO TIME (ns)
# print(data.groupby(['CHAIN']).mean())
# pivot_chainA = pd.pivot_table(data[data['CHAIN']=='A'],columns='ID',index='TIME',values='COUNT')
# # print(pivot)
# pivot_chainB = pd.pivot_table(data[data['CHAIN']=='B'],columns='ID',index='TIME',values='COUNT')

# fig,axes = plt.subplots(1,2)

# plot1 = sns.heatmap(pivot_chainA,ax=axes[0],cmap='viridis',vmin=0,vmax=5,cbar_kws={'label':"Number of PIP2 Binding"})
# plot2 = sns.heatmap(pivot_chainB,ax=axes[1],cmap='viridis',vmin=0,vmax=5,cbar_kws={'label':"Number of PIP2 Binding"})

# plot1.set_xlabel("Residue")
# plot2.set_xlabel("Residue")
# plot1.set_ylabel("Time (ns)")
# plot2.set_ylabel("Time (ns)")

# plot1.locator_params(axis='y', nbins=8)
# plot2.locator_params(axis='y', nbins=8)

# # # ### MISCELLANEOUS ###
# plt.suptitle("%s"%(ifile[:-4]),va='top')
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,5)   ## Wide x Height
# plt.tight_layout()
# # plt.savefig("%s"%(ifile[:-3]))
# plt.show()
