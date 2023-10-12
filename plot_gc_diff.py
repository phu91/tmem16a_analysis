import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc

sns.set_context("talk")

### FUNCTION 
def read_gc_output(input,nres):
    data = []
    with open(input,"r+") as ifile:
        lines = ifile.readlines()
        for line in lines:
            line = line.split()
            for each in line:
                data.append(each)
    new_data = np.reshape(data, (nres, nres)).T
    return new_data.astype(float)

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--i1', type=str, default='',
                    help='INPUT to GC 1')

parser.add_argument('--i2', type=str, default='',
                    help='INPUT to GC 2')

parser.add_argument('--nr', type=int, default='',
                    help='Number of Residues')

parser.add_argument('--sn', type=str, default='UNKNOWN',
                    help='System Name')

args = parser.parse_args()

ifile1 =  args.i1
ifile2 =  args.i2
nres   =  args.nr
sysname = args.sn

chain1 = read_gc_output(ifile1,nres)
chain2 = read_gc_output(ifile2,nres)

diff = abs(chain1-chain2)
##########
mask = np.triu(diff)

cmap = sns.diverging_palette(150,20, as_cmap=True)

fig, ax = plt.subplots(1)
ax = sns.heatmap(diff,
            cmap='icefire',
            xticklabels=50,
            yticklabels=50,
            # robust=True,
            # mask=mask,
            cbar=True,
            cbar_kws={"label":"GC(A)-GC(B)"},
            square=True,
            center=0,
            ax=ax,
)
# ax.tick_params(top=True,labeltop=True,bottom=False,labelbottom=False)
cbar = ax.collections[0].colorbar
cbar.set_label('Delta(GC)', rotation=-90,labelpad=20)

plt.xticks(rotation=45)
plt.ylim(0,nres)
plt.xlim(0,nres)

# # ### MISCELLANEOUS ###
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(10,7.5)
# plt.locator_params(axis='both', nbins=798)
plt.tight_layout()
plt.savefig("%s"%("GCorrelation-%s.png"%(sysname)))
plt.show()
