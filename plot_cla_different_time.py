import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
import colorcet as cc
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker

sns.set(style="ticks", context="notebook")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT Chloride profile')

parser.add_argument('--topol', type=str, default='',
                    help='INPUT PORE Topology')

parser.add_argument('--factor', type=float, default='1',
                    help='Conversion factor from FRAME to NS. Default 1 frame/ns.'
                    )

args = parser.parse_args()

ifile =  args.input
pore_top = args.topol
factor = args.factor

data = pd.read_csv(ifile,comment='#',
delim_whitespace=True,
names=['FRAME','INDEX','CHAIN','ZPOS','POSITION_PORE','DCOM_PORE'])

data['ns']=data.FRAME*factor

topol_data = pd.read_csv(pore_top,comment='#',
delim_whitespace=True,
names=['PART','CHAIN','MIN_Z','MAX_Z'])
topol_data=topol_data.groupby(['PART']).mean().reset_index()

# print(topol_data)
### DEFINE PORE TOPOLOGY
Ca_Binding=topol_data.loc[topol_data.PART=="Ca_binding"]
Inner_Vestibule=topol_data.loc[topol_data.PART=="Inner_Vest"]
Neck=topol_data.loc[topol_data.PART=="Neck"]
Outer_Vestibule=topol_data.loc[topol_data.PART=="Outer_Vest"]

# print(Outer_Vestibule)
pore_topology_color = sns.color_palette("Accent_r",4)
pore_topology_label = ['Ca2+ Binding','Inner Vestibule','Neck','Outer Vestibule']
pore_topology_name = [Ca_Binding,Inner_Vestibule,Neck,Outer_Vestibule]

fig, axes = plt.subplots(1,2,sharey=True)
marker_colors_a = ListedColormap(cc.glasbey_warm[:20])
marker_colors_b = ListedColormap(cc.glasbey_cool[:20])
color_list=[marker_colors_a,marker_colors_b]
chain_list=['A','B']
marker_list=['o','X']
weight='bold'
fontsize=15

min_z_plot = np.round(topol_data.MIN_Z.min()-30,0)
max_z_plot = np.round(topol_data.MAX_Z.max()+30,0)

if len(data.ns)==0:
    for ind,ax in enumerate(axes):
        for i,j,k in zip(pore_topology_name,pore_topology_label,pore_topology_color):
            ax.fill_between([0,100000],i.MIN_Z,i.MAX_Z,color=k,label=j,alpha=0.6,lw=0)
        ax.set_xticklabels(())
        ax.set_xlabel("Time (ns)",weight=weight,fontsize=fontsize)
        ax.set_ylabel("Z Direction",weight=weight,fontsize=fontsize)
        ax.set_ylim([min_z_plot,max_z_plot])
        ax.set_title("CHAIN %s"%(chain_list[ind]),weight=weight,fontsize=fontsize)

else:
    for ind,ax in enumerate(axes):
        for i,j,k in zip(pore_topology_name,pore_topology_label,pore_topology_color):
            ax.fill_between([0,5000],i.MIN_Z,i.MAX_Z,color=k,label=j,alpha=0.6,lw=0)
        g = sns.scatterplot(data=data.query("CHAIN=='%s'"%(chain_list[ind])),
                    x="ns",
                    y="ZPOS",
                    hue="INDEX",
                    # style="CHAIN",
                    alpha=0.85,
                    marker=marker_list[ind],
                    # style_order=['A','B'],
                    legend=None,
                    palette=color_list[ind],
                    linewidth=0.8,
                    edgecolor='gray',
                    s=50,
                    ax=ax
                    )
        ax.set_xlabel("Time (ns)",weight=weight,fontsize=fontsize)
        ax.set_ylabel("Z Direction",weight=weight,fontsize=fontsize)
        ax.set_ylim([min_z_plot,max_z_plot])
        ax.xaxis.set_major_locator(ticker.AutoLocator())
        # ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        ax.set_title("CHAIN %s"%(chain_list[ind]),weight=weight,fontsize=fontsize)
# # # ### MISCELLANEOUS ###
# plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,5.5)
plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()

plt.savefig("%s.png"%(ifile[:-4]),dpi=700)

# plt.show()
