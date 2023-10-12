# Updated on 11/29/2022
import pandas as pd
import numpy as np
import os 
import warnings
warnings.filterwarnings("ignore")
import re
import sys
import matplotlib.pyplot as plt
import seaborn as sns 
from mpl_toolkits.axes_grid1 import make_axes_locatable

sns.set_context("paper")

histo_1 = sys.argv[1]
histo_2 = sys.argv[2]
cla_pos = sys.argv[3]
pore_def = sys.argv[4]

df1 = pd.read_csv(histo_1,
                  comment='#',
                  delim_whitespace=True,
                  names=['edge1','edge2','half_edge','mean','std'])

df2 = pd.read_csv(histo_2,
                  comment='#',
                  delim_whitespace=True,
                  names=['edge1','edge2','half_edge','mean','std'])

df3 = pd.read_csv(cla_pos,
                  comment='#',
                  delim_whitespace=True,
                  names=['f','index','z','chain'])

df4 = pd.read_csv(pore_def,
                  comment='#',
                  delim_whitespace=True,
                  names=['f','ca_min','ca_max','inner_min','inner_max','neck_min','neck_max','outer_min','outer_max'])

df5 = df4.mean().drop(['f'])
# print(df5)
##################
lw = 3

fig = plt.figure()
ax = fig.add_subplot(111)
sns.lineplot(data=df1,
             x='half_edge',
             y='mean',
             ax=ax,
             label="A",
             linewidth=lw
             )

std_up = df1['mean']+df1['std']
std_down = df1['mean']-df1['std']
ax.fill_between(df1['half_edge'],std_up,std_down,alpha=0.2)

sns.lineplot(data=df2,
             x='half_edge',
             y='mean',
             ax=ax,
             label="B",
             linewidth=lw
             )
std_up = df2['mean']+df2['std']
std_down = df2['mean']-df2['std']
ax.fill_between(df2['half_edge'],std_up,std_down,alpha=0.2)

cl_a_list=['148898']
for i in np.arange(len(cl_a_list)):
	# print(cl_b_list[i])
    sns.scatterplot(data=df3.query("index==%s"%(cl_a_list[i])),
                    x='z',
                    y=-6-(i*2),
                    size='f',
                    sizes=(10,300),
                    hue='f',
                    # style='index',
                    marker='^',
                    legend=None,
                    ax=ax,
                    # alpha=0.8,
                    palette='Blues')

cl_b_list=['149256','148998']
for i in np.arange(len(cl_b_list)):
	# print(cl_b_list[i])
    sns.scatterplot(data=df3.query("index==%s"%(cl_b_list[i])),
                    x='z',
                    y=-10-(i*2),
                    size='f',
                    sizes=(10,300),
                    hue='f',
                    # style='index',
                    marker='^',
                    legend=None,
                    ax=ax,
                    # alpha=0.8,
                    palette='Oranges')

ax.scatter(df3.z[0],df3.z[1],c='black',s=200,label='Chloride',marker='^')

timeline = np.arange(0,2,1)
sm = plt.cm.ScalarMappable(cmap="binary", norm=None)
sm.set_array([timeline.min(),timeline.max()])
ax.get_legend().remove()
cbar = ax.figure.colorbar(sm,location='top',orientation='horizontal', pad=0.01,aspect=90,ticks=[0, 1])
cbar.ax.set_xticklabels(['Start', 'End'])
cbar.set_label("Color Intensity shows the Timepoint of Chlorides",
               rotation=0,
               labelpad=1,
               fontsize=10)
               
# sm = plt.cm.ScalarMappable(cmap="Orange", norm=None)
# sm.set_array([df3.f.min(), df3.f.max()])
# ax.get_legend().remove()
# cbar = ax.figure.colorbar(sm,location='top',orientation='horizontal', pad=0.01,aspect=90)
# cbar.set_label("Frame",
#                rotation=0,
#                labelpad=1,
#                fontsize=10)

a = list(df5.values)
pore_label = ["Ca2+ Binding","Inner Vestibule","Neck","Outer Vestibule"]
color_pore = sns.color_palette('Set2')
offset = np.linspace(10,20,4)

for i,c in zip(np.arange(len(pore_label)),color_pore):
    ax.plot((a[2*i],a[2*i+1]),(offset[i],offset[i]),color=c,linewidth=lw,label=pore_label[i],marker='s')


plt.locator_params(axis='y', nbins=10)
plt.xlim([-21,170])
plt.xlabel("Z Coordinate",fontweight='bold',size=20)
plt.ylabel("Average Radii",fontweight='bold',size=20)
plt.legend(loc='upper right',
        #    , bbox_to_anchor=(0.5, 1.05),
        #   ncol=3, fancybox=True, shadow=True
        )

plt.tight_layout()
plt.savefig("pore_profile_sys2_rep2.png",dpi=300)
plt.show()