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

chaina = read_gc_output(ifile1,nres)
chainb = read_gc_output(ifile2,nres)

# print(chaina)
##########

fig, ax = plt.subplots(1)

maskhigh = np.triu(chaina)
g1 = sns.heatmap(chaina,
            cmap="YlGn",
            xticklabels=50,
            yticklabels=50,
            robust=True,
            mask=maskhigh,
            cbar_kws={"orientation": "horizontal","location":"top","pad": 0.02,"fraction":0.036,"label":"CHAIN A"},
            ax=ax,
            )

masklow = np.tril(chainb)
g2 = sns.heatmap(chainb,
            cmap="PuRd",
            xticklabels=50,
            yticklabels=50,
            robust=True,
            mask=masklow,
            cbar_kws={"orientation": "vertical","location":"right","pad": 0.02,"label":"CHAIN B"},
            ax=ax,
            )

plt.xticks(rotation=45) 
plt.ylim(0,nres)
plt.xlim(0,nres)

# ### MISCELLANEOUS ###
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(8,7.5)
# plt.locator_params(axis='both', nbins=798)
plt.tight_layout()
plt.savefig("%s"%("GCorrelation-%s.png"%(sysname)))
plt.show()
