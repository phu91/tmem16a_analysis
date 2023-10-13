# Updated on 11/29/2022
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math, sys
import seaborn as sns 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
import MDAnalysis as mda
from matplotlib import rc

sns.set(style="ticks", context="talk")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to PIP2 positions')

parser.add_argument('--helix', type=str, default='',
                    help='Last frame of PRODUCTION to extract the helix positions')

args = parser.parse_args()

ifile =  args.input
helix_pos = args.helix

u = mda.Universe(helix_pos,helix_pos)

chainA_com = u.select_atoms("segid PROA",updating=True).center_of_geometry()
chainB_com = u.select_atoms("segid PROB",updating=True).center_of_geometry()

helix_01a = u.select_atoms("segid PROA and resid 327 to 360").center_of_geometry()
helix_02a = u.select_atoms("segid PROA and resid 399 to 439").center_of_geometry()
helix_03a = u.select_atoms("segid PROA and resid 493 to 518").center_of_geometry()
helix_04a = u.select_atoms("segid PROA and resid 534 to 562").center_of_geometry()
helix_05a = u.select_atoms("segid PROA and resid 572 to 598").center_of_geometry()
helix_06a = u.select_atoms("segid PROA and resid 629 to 664").center_of_geometry()
helix_07a = u.select_atoms("segid PROA and resid 692 to 713").center_of_geometry()
helix_08a = u.select_atoms("segid PROA and resid 714 to 738").center_of_geometry()
helix_09a = u.select_atoms("segid PROA and resid 753 to 775").center_of_geometry()
helix_10a = u.select_atoms("segid PROA and resid 854 to 881").center_of_geometry()

helix_01b = u.select_atoms("segid PROB and resid 327 to 360").center_of_geometry()
helix_02b = u.select_atoms("segid PROB and resid 399 to 439").center_of_geometry()
helix_03b = u.select_atoms("segid PROB and resid 493 to 518").center_of_geometry()
helix_04b = u.select_atoms("segid PROB and resid 534 to 562").center_of_geometry()
helix_05b = u.select_atoms("segid PROB and resid 572 to 598").center_of_geometry()
helix_06b = u.select_atoms("segid PROB and resid 629 to 664").center_of_geometry()
helix_07b = u.select_atoms("segid PROB and resid 692 to 713").center_of_geometry()
helix_08b = u.select_atoms("segid PROB and resid 714 to 738").center_of_geometry()
helix_09b = u.select_atoms("segid PROB and resid 753 to 775").center_of_geometry()
helix_10b = u.select_atoms("segid PROB and resid 854 to 881").center_of_geometry()

pip2 = u.select_atoms("resname PLPI24 and name P",updating=True)

df = pd.read_csv(ifile,
                 comment='#',
                 names=['f','id','x','y','d_a'],
                 delim_whitespace=True)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
markersize=20

text_params = {'ha': 'right','fontweight': 'bold','color': 'white','fontsize':'20',}

ax.plot(helix_01a[0],helix_01a[1],marker="d",markersize=markersize,color='lightpink',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_02a[0],helix_02a[1],marker="d",markersize=markersize,color='lightpink',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_03a[0],helix_03a[1],marker="d",markersize=markersize,color='lightpink',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_04a[0],helix_04a[1],marker="d",markersize=markersize,color='lightpink',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_05a[0],helix_05a[1],marker="d",markersize=markersize,color='lightpink',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_06a[0],helix_06a[1],marker="d",markersize=markersize,color='lightpink',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_07a[0],helix_07a[1],marker="d",markersize=markersize,color='lightpink',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_08a[0],helix_08a[1],marker="d",markersize=markersize,color='lightpink',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_09a[0],helix_09a[1],marker="d",markersize=markersize,color='lightpink',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_10a[0],helix_10a[1],marker="d",markersize=markersize,color='lightpink',markeredgewidth=1, markeredgecolor='white')

ax.text(helix_01a[0],helix_01a[1],"1A",**text_params)
ax.text(helix_02a[0],helix_02a[1],"2A",**text_params)
ax.text(helix_03a[0],helix_03a[1],"3A",**text_params)
ax.text(helix_04a[0],helix_04a[1],"4A",**text_params)
ax.text(helix_05a[0],helix_05a[1],"5A",**text_params)
ax.text(helix_06a[0],helix_06a[1],"6A",**text_params)
ax.text(helix_07a[0],helix_07a[1],"7A",**text_params)
ax.text(helix_08a[0],helix_08a[1],"8A",**text_params)
ax.text(helix_09a[0],helix_09a[1],"9A",**text_params)
ax.text(helix_10a[0],helix_10a[1],"10A",**text_params)

ax.plot(helix_01b[0],helix_01b[1],marker="o",markersize=markersize,color='skyblue',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_02b[0],helix_02b[1],marker="o",markersize=markersize,color='skyblue',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_03b[0],helix_03b[1],marker="o",markersize=markersize,color='skyblue',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_04b[0],helix_04b[1],marker="o",markersize=markersize,color='skyblue',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_05b[0],helix_05b[1],marker="o",markersize=markersize,color='skyblue',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_06b[0],helix_06b[1],marker="o",markersize=markersize,color='skyblue',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_07b[0],helix_07b[1],marker="o",markersize=markersize,color='skyblue',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_08b[0],helix_08b[1],marker="o",markersize=markersize,color='skyblue',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_09b[0],helix_09b[1],marker="o",markersize=markersize,color='skyblue',markeredgewidth=1, markeredgecolor='white')
ax.plot(helix_10b[0],helix_10b[1],marker="o",markersize=markersize,color='skyblue',markeredgewidth=1, markeredgecolor='white')

ax.text(helix_01b[0],helix_01b[1],"1B",**text_params)
ax.text(helix_02b[0],helix_02b[1],"2B",**text_params)
ax.text(helix_03b[0],helix_03b[1],"3B",**text_params)
ax.text(helix_04b[0],helix_04b[1],"4B",**text_params)
ax.text(helix_05b[0],helix_05b[1],"5B",**text_params)
ax.text(helix_06b[0],helix_06b[1],"6B",**text_params)
ax.text(helix_07b[0],helix_07b[1],"7B",**text_params)
ax.text(helix_08b[0],helix_08b[1],"8B",**text_params)
ax.text(helix_09b[0],helix_09b[1],"9B",**text_params)
ax.text(helix_10b[0],helix_10b[1],"10B",**text_params)

g = sns.kdeplot(
    data=df.query("-60<x<60 and -60<y<60"),
    x="x",
    y="y",
    # hue="kind",
    cmap='plasma', 
    shade=True,
    # cbar=True,
    thresh=0,
    #levels=100
)

# cax = fig.add_axes([ax.get_position().x1+0.05,ax.get_position().y0,0.02,ax.get_position().height])
# cb = fig.colorbar(hb,cax=cax)
# cb.set_label('Frames')

# ax.set_xlim([20,150])
# ax.set_ylim([20,150])

ax.set_xlabel('X')
ax.set_ylabel('Y')

# ### MISCELLANEOUS ###
plt.xlim([-60,60])
plt.ylim([-60,60])
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(10,7.5)
plt.locator_params(axis='y', nbins=6)
# plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()
plt.savefig("%s"%(ifile[:-3]))
plt.show()
