import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc

sns.set_context("notebook")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to RMSF profile')

args = parser.parse_args()

ifile =  args.input

data = pd.read_csv(ifile,
                 names=['Frame','Chain','Helix','pip2'],
                 comment='#',
                 delim_whitespace=True)

# print(data)

chain_list = ['A','B']
color_list = sns.color_palette('Set2')

data2 = data.groupby(['Frame','Chain']).sum()
data2.to_csv("tmp")
df1 = pd.read_csv("tmp")

fig,axes = plt.subplots(2,1,sharex=False,sharey=False,layout="constrained")
axes = axes.flatten()

data['pip2']=data['pip2']/60*100

g1 = sns.boxplot(data=data,
            x='Helix',
            y='pip2',
            ax=axes[0],
            hue='Chain',
            # split=True,
            )

g1.set_ylabel("% PIP2")

df1['pip2']=df1['pip2']/60*100

g2 = sns.boxplot(data=df1,
            x='Chain',
            y='pip2',
            ax=axes[1]
            # hue='Chain',
            # split=True,
            )
g2.set_ylabel("% PIP2")


# ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,10)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
plt.savefig("%s"%(ifile[:-3]))
plt.show()
