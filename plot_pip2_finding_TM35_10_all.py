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

parser.add_argument('--filter', nargs='+', help='<Required> Set flag', required=True)

args = parser.parse_args()

ifile =  args.input
filter= args.filter

data = pd.read_csv(ifile,comment='#',
delim_whitespace=True,names=['frame','resid','position','binding']
)

data['time']=data['frame']*0.2
### COLECTING INFORMATION FROM FILTER
residue_count = len(filter)//2
residue_id = []
residue_position = []
for i in range(len(filter)//2):
    # print(filter[i+i*i])
    residue_id.append(filter[i*2])
    residue_position.append(filter[i*2+1])

filtered_df = pd.DataFrame()
for i in range(residue_count):
    filtered = data.query("resid==%s and position=='%s'"%(residue_id[i],residue_position[i]))
    filtered.loc[filtered['position']=='567B_887A','position']="TM10_A"
    filtered.loc[filtered['position']=='567A_887B','position']="TM10_B"
    filtered_df = pd.concat([filtered,filtered_df],ignore_index=True)
# print(filtered_df)
filtered_df = filtered_df.sort_values(by=['position'])
# fig,axes = plt.subplots(1,2,sharex=True,sharey=True)
# marker_colors = cc.glasbey[:2]
# print(marker_colors)

position_name_list=['TM10_A','TM10_B']
g = sns.violinplot(data=filtered_df,
x='position',
y='binding',)
g.set_yticks(np.arange(0, 1.1, 1))
g.set_yticklabels(['No','Yes'])


# for ind,ax in enumerate(axes.flatten()):
# #     # print(ind,ax)
# #     # print(residue_id[ind],residue_position[ind])
#     g = sns.boxplot(data=filtered_df.query("position=='%s'"%(position_name_list[ind])),
#     x='time',
#     y='binding',
#     # drawstyle='steps-post',
# #     # marker='o',
#     color=marker_colors[ind],
#     ax=ax
#     )
#     g.set_ylim([-0.1,1.1])
#     g.set_title("PIP2 || BINDING %s"%(position_name_list[ind]))
#     g.set_ylabel("Binding")
#     g.set_xlabel("Time(ns)")
#     g.set_yticks(np.arange(0, 1.1, 1))
#     g.set_yticklabels(['No','Yes'])

plt.suptitle("%s"%(ifile))
# # plt.rcParams['ps.useafm'] = True
# # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# # plt.rcParams['pdf.fonttype'] = 42
# # plt.gcf().set_size_inches(7.5,5.5)
# # plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()
# # plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
plt.show()