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

data = pd.read_csv(ifile,comment='#',
                   delim_whitespace=True,
                   names=['Frame','PC1','PC2','PC3','Chain','ps'])

data['ps'] = data['ps']/1000

data = data.astype({'ps':'int'})
# print(data['ps'][::1000])

chain_list = ['A','B','C']
color_list = sns.color_palette('Set2')


fig,axes = plt.subplots(3,1,sharex=True,sharey=True,layout="constrained")
axes = axes.flatten()

for ax in enumerate(axes):
    # print(data.query("Chain=='%s'"%(chain_list[ax[0]])))
    g = sns.kdeplot(data=data.query("Chain=='%s'"%(chain_list[ax[0]])),
            x='PC1',
            y='PC2',
            fill=True,
            # alpha=0.75,
            color=color_list[ax[0]],
            ax=ax[1]
            )
    g = sns.scatterplot(data=data.query("Chain=='%s'"%(chain_list[ax[0]]))[::1000],
            x='PC1',
            y='PC2',
            hue='ps',
            ax=ax[1],
                  )
    g.legend(title="Chain %s"%(chain_list[ax[0]]))


    
# ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,6)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
plt.savefig("%s"%(ifile[:-3]))
plt.show()
