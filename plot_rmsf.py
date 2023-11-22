import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc

<<<<<<< HEAD
<<<<<<< HEAD
sns.set_context("talk")
=======
sns.set_context("notebook")
>>>>>>> pip2
=======
sns.set_context("talk")
>>>>>>> 412b5a255f3017dd65aa1dd65a64a731574bb88b

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to RMSF profile')

args = parser.parse_args()

ifile =  args.input

data = pd.read_csv(ifile,comment='#',
                   delim_whitespace=True,
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 412b5a255f3017dd65aa1dd65a64a731574bb88b
                   names=['helix','chain','resid','rmsf'])


fig,axes = plt.subplots(2,5,sharex=True,sharey=True,squeeze=True)
axes = axes.flatten()

for ax in enumerate(axes):
    # print(ax)
    sns.lineplot(x='resid',y='rmsf',data=data.query("helix==%s"%(ax[0]+1)),
                 hue='chain',
                 palette='Set2',
                 ax=ax[1])
    ax[1].set_title("TM %s"%(ax[0]+1))

# ### MISCELLANEOUS ###
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(18,7.5)
plt.locator_params(axis='both', nbins=5)
<<<<<<< HEAD
=======
                   names=['chain','resname','resid','rmsf','system'])


colors = sns.color_palette("husl", 10)

system_list = ['0A0B_R1','3A3B_R1','3A0B_R1','3A1B_R1','0A0B_R2','3A3B_R2','3A0B_R2','3A1B_R2']
# print(data)
fig,axes = plt.subplots(2,4,sharex=True,sharey=True)
axes = axes.flatten()

for ind,ax in enumerate(axes):
    # print(ind,ax)
    sns.lineplot(data=data.query("system=='%s'"%(system_list[ind])), 
                x='resid', y='rmsf',
                hue='chain',
                palette='Set2',
                lw=3,
                ax=ax)
    ax.axvspan(327, 360, zorder=0, alpha=0.2, label='TM 01',color=colors[0])
    ax.axvspan(399, 439, zorder=0, alpha=0.2, label='TM 02',color=colors[1])
    ax.axvspan(478, 518, zorder=0, alpha=0.2, label='TM 03',color=colors[2])
    ax.axvspan(534, 562, zorder=0, alpha=0.2, label='TM 04',color=colors[3])
    ax.axvspan(572, 598, zorder=0, alpha=0.2, label='TM 05',color=colors[4])
    ax.axvspan(629, 664, zorder=0, alpha=0.2, label='TM 06',color=colors[5])
    ax.axvspan(692, 713, zorder=0, alpha=0.2, label='TM 07',color=colors[6])
    ax.axvspan(718, 740, zorder=0, alpha=0.2, label='TM 08',color=colors[7])
    ax.axvspan(753, 780, zorder=0, alpha=0.2, label='TM 09',color=colors[8])
    ax.axvspan(854, 927, zorder=0, alpha=0.2, label='TM 10',color=colors[9])
    ax.set_title("%s"%(system_list[ind]))

# ### MISCELLANEOUS ###
# plt.legend()
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(14,8)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
>>>>>>> pip2
=======
>>>>>>> 412b5a255f3017dd65aa1dd65a64a731574bb88b
plt.tight_layout()
plt.savefig("%s"%(ifile[:-3]))
plt.show()
