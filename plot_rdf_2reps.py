import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc

sns.set_context("talk")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to RDF profile')

args = parser.parse_args()

ifile =  args.input

data = pd.read_csv(ifile,comment='#',
                   delim_whitespace=True,
                   names=['rep','helix','chain','radius','gofR'])

# print(data.head)
# palette = sns.color_palette("icefire", as_cmap=True)

fig,axes = plt.subplots(2,5,sharex=True,squeeze=True)
axes = axes.flatten()

for ax in enumerate(axes):
    # print(ax)
    sns.lineplot(x='radius',y='gofR',data=data.query("helix==%s"%(ax[0]+1)),
                 hue='chain',
                #  style='rep',
                 palette='Set2',
                 ax=ax[1])
    ax[1].set_title("TM %s"%(ax[0]+1))


# ### MISCELLANEOUS ###
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(18,7.5)
plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
plt.savefig("%s"%(ifile[:-3]))
plt.show()
