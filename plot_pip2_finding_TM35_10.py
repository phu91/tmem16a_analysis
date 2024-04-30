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

# print(residue_id)
# print(residue_position)
fig,axes = plt.subplots(residue_count,1,sharex=True,sharey=True)
marker_colors = cc.glasbey[:int(residue_count)]

for ind,ax in enumerate(axes.flatten()):
    # print(ind,ax)
    # print(residue_id[ind],residue_position[ind])
    g = sns.lineplot(data=data.query("resid==%s and position=='%s'"%(residue_id[ind],residue_position[ind])),
    x='time',
    y='binding',
    drawstyle='steps-post',
    # marker='o',
    color=marker_colors[ind],
    ax=ax
    )
    g.set_ylim([-0.1,1.1])
    g.set_title("PIP2 || RESID %s || BINDING %s"%(residue_id[ind],residue_position[ind]))
    g.set_ylabel("Binding")
    g.set_xlabel("Time(ns)")
    g.set_yticks(np.arange(0, 1.1, 1))
    g.set_yticklabels(['No','Yes'])

# plt.suptitle("%s"%(ifile))
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,5.5)
# plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()
# plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
plt.show()