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
                   names=['FRAME','CHAIN','RESIDUE','RESID','OMEGA','PHI','PSI','CHI1','CHI2','SYSTEM'])

res_list=['879','884']
dihedral_list=['CHI1','CHI2']
# fig,axes = plt.subplots(3,2,sharex=True,sharey=False,layout="constrained")
# axes = axes.flatten()

# for row in range(len(res_list)):
#     for col in range(len(dihedral_list)):
#         g = sns.scatterplot(data=data.query("RESID==%s & DIHEDRAL=='%s'"%(res_list[row],dihedral_list[col])),
#                     x='FRAME',
#                     y='ANGLE',
#                     hue='CHAIN',
#                     style='CHAIN',
#                     ax=axes[row][col],
#                     alpha=0.5,
#                     lw=0,
#                     )
#         g.set_ylabel("%s RESID %s"%(dihedral_list[col],res_list[row]))

#############################
fig,axes = plt.subplots(1,2,sharex=True,sharey=False,layout="constrained")
axes = axes.flatten()
cmap = sns.color_palette('Set2')
for col in range(len(res_list)):
    g = sns.kdeplot(data=data.query("RESID==%s & CHAIN=='B'"%(res_list[col])),
                x='CHI1',
                y='CHI2',
                # hue='CHAIN',
                fill=True,
                alpha=0.90,
                color=cmap[1],
                # palette='Set2',
                # cbar=True,
                # z=1,
                ax=axes[col],
                )
    g = sns.kdeplot(data=data.query("RESID==%s & CHAIN=='A'"%(res_list[col])),
                x='CHI1',
                y='CHI2',
                # hue='CHAIN',
                fill=True,
                alpha=0.90,
                color=cmap[0],
                # palette='Set2',
                # z=0,
                # cbar=True,
                ax=axes[col],
                )
    # g = sns.rugplot(data=data.query("RESID==%s"%(res_list[col])),
    #             x='CHI1',
    #             y='CHI2',
    #             hue='CHAIN',
    #             palette='Set2',
    #             ax=axes[col],
    #             )
    g.set_title("ASP %s"%(res_list[col]))
# color_list = sns.color_palette('Set2')


# ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,4)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
# plt.tight_layout()
plt.savefig("KDE%s"%(ifile[:-3]))
plt.show()
