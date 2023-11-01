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
                   names=['FRAME','CHAIN','RESIDUE NAME','RESID','OMEGA','PHI','PSI','SYSTEM'])
                   
# print(data)
sns.lineplot(data=data,
            x='FRAME',
            y='OMEGA',
            style='CHAIN',
            hue='RESID',
            # markers=True,
            palette='Set1',
            )

color_list = sns.color_palette('Set2')


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
