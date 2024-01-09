import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
import colorcet as cc

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

### DEFINE PORE TOPOLOGY
PORE = data[['POSITION_PORE','ZPOS']]
Inner_Vestibule=PORE.loc[PORE.POSITION_PORE=="Inner_Vestibule"]
Ca_Binding=PORE.loc[PORE.POSITION_PORE=="Ca_Binding"]
Neck=PORE.loc[PORE.POSITION_PORE=="Neck"]
Outer_Vestibule=PORE.loc[PORE.POSITION_PORE=="Outer_Vestibule"]
Out_of_Zone=PORE.loc[PORE.POSITION_PORE=="Out_of_Zone"]

# print(Outer_Vestibule.min()[1])
pore_topology_color = sns.color_palette("Set3")
pore_topology_label = ['Inner Vestibule','Ca2+ Binding','Neck','Outer Vestibule']
pore_topology_name = [Inner_Vestibule,Ca_Binding,Neck,Outer_Vestibule]

fig, ax = plt.subplots(1,1)
marker_colors = sns.color_palette(cc.glasbey_bw, n_colors=20)

for i,j,k in zip(pore_topology_name,pore_topology_label,pore_topology_color):
    ax.fill_between(data.query("DCOM_PORE<=%s"%(cutoff)).ns,i.min()[1],i.max()[1],color=k,label=j)
g = sns.scatterplot(data=data.query("DCOM_PORE<=%s"%(cutoff)), 
            x="ns", 
            y="ZPOS", 
            hue="INDEX",
            style="CHAIN",
            alpha=0.85,
            markers=True,
            style_order=['A','B'],
            legend=None,
            palette=marker_colors,
            linewidth=0.8,
            edgecolor='white',
            s=80,
            ax=ax
            )
weight='bold'
fontsize=20

ax.set_xlabel("Time (ns)",weight=weight,fontsize=fontsize)
ax.set_ylabel("Z Coordinate",weight=weight,fontsize=fontsize)
ax.set_ylim([Ca_Binding.min()[1],Outer_Vestibule.max()[1]])
ax.xaxis.set_tick_params(labelsize=fontsize-5)
ax.yaxis.set_tick_params(labelsize=fontsize-5)
# ax.set_xlim([0,1000])

# # ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,5.5)
plt.legend(loc='upper left',bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()
plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
# plt.show()
