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
data_plot['time']=data_plot['frame']*0.2

data_plot_pivot = data_plot.pivot_table(index='time',columns='position',values='binding')
data_plot_pivot = data_plot_pivot.reset_index()
# print(data_plot_pivot)
data_plot_pivot = data_plot_pivot[['time','TM10_A','TM10_B']]
# print(data_plot_pivot)
data_plot_pivot.to_csv("%s_PLOT.dat"%(ifile[:-4]),sep='\t',index=False)

fig, axes = plt.subplots(2,1)

g1 = sns.violinplot(data=data_plot,
x='position',
y='binding',
ax=axes[0])
g1.set_ylabel("Binding Mode")
g1.set_xlabel("Binding Position")
g1.set_yticks(np.arange(0, 1.5, 1))
g1.set_yticklabels(['No','Yes'])

g2 = axes[1].plot(data_plot_pivot.time,data_plot_pivot.TM10_A,data_plot_pivot.time,data_plot_pivot.TM10_B,)
axes[1].set_xlabel("Time(ns)")
axes[1].set_ylabel("Binding Mode")
axes[1].set_yticks(np.arange(0, 1.5, 1))
axes[1].set_yticklabels(['No','Yes'])


plt.suptitle("%s"%(ifile))
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,5.5)
plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()
# plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
plt.show()