import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
import colorcet as cc

sns.set(style="ticks", context="notebook")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT Chloride profile')

args = parser.parse_args()

ifile =  args.input

data = pd.read_csv(ifile,comment='#',
delim_whitespace=True,names=['frame','position','binding']
)

# data['time']=data['frame']*0.2
# print(data)
position_count = len(pd.unique(data['position']))
# print(position_count)
position_current = pd.unique(data['position'])
# print(position_current)
with open (ifile,'a') as fileNow:
    if position_count<2:
        if position_current=='TM10_A':
            fileNow.write("0    000 TM10_B  0\n")
        else:
            fileNow.write("0    000 TM10_A  0\n")

data_plot = pd.read_csv(ifile,comment='#',
delim_whitespace=True,names=['frame','position','binding']
)
# print(data_plot)
data_plot_pivot = data_plot.pivot_table(index='frame',columns='position',values='binding')
data_plot_pivot = data_plot_pivot.reset_index()
# print(data_plot_pivot)
data_plot_pivot = data_plot_pivot[['frame','TM10_A','TM10_B']]
# print(data_plot_pivot)
data_plot_pivot.to_csv("%s_PLOT.dat"%(ifile[:-4]),sep='\t',index=False)
# print(data_plot_pivot.columns)

# fig,axes = plt.subplots(1,2,sharex=True,sharey=True)
# marker_colors = cc.glasbey[:2]
# # print(marker_colors)

# g = sns.violinplot(data=data_plot_pivot,
# x='position',
# y='binding',)

# # g.set_ylim([-0.1,1.1])
# g.set_ylabel("Binding Mode")
# g.set_xlabel("Binding Position")
# g.set_yticks(np.arange(0, 1.1, 1))
# g.set_yticklabels(['No','Yes'])

# plt.suptitle("%s"%(ifile))
# # plt.rcParams['ps.useafm'] = True
# # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# # plt.rcParams['pdf.fonttype'] = 42
# # plt.gcf().set_size_inches(7.5,5.5)
# # plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
# plt.tight_layout()
# # plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
plt.show()