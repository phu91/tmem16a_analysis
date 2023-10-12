import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc

sns.set_context("talk")

sys1 = pd.read_csv("RMSF_SEGid_withCa-NoPIP2.dat",comment='#',
                   delim_whitespace=True,
                   names=['chain','resid','rmsf','system'])

sys2 = pd.read_csv("RMSF_SEGid_withCa-withPIP2_rep1.dat",comment='#',
                   delim_whitespace=True,
                   names=['chain','resid','rmsf','system'])

sys3 = pd.read_csv("RMSF_SEGid_noCa-noPIP2.dat",comment='#',
                   delim_whitespace=True,
                   names=['chain','resid','rmsf','system'])

sys4 = pd.read_csv("RMSF_SEGid_noCa-withPIP2_rep1.dat",comment='#',
                   delim_whitespace=True,
                   names=['chain','resid','rmsf','system'])

data = [sys1,sys2,sys3,sys4]

fig,axes = plt.subplots(4,1,
                        sharex=True,
                        sharey=True,
                        squeeze=True,
                        gridspec_kw = {'wspace':0.2, 'hspace':0.2})

point_color = sns.cubehelix_palette(10,start=.5, rot=20, as_cmap=True)

for ax in enumerate(axes):
    # print(ax)
    sns.lineplot(x='resid',y='rmsf',data=data[ax[0]],
                 hue='chain',
                 palette='Set2',
                 ax=ax[1],
                 legend=False,
                 )

    ax[1].axvspan(327, 360, zorder=0, alpha=0.2, label='TM 01')
    ax[1].axvspan(399, 439, zorder=0, alpha=0.2, label='TM 02')
    ax[1].axvspan(493, 518, zorder=0, alpha=0.2, label='TM 03')
    ax[1].axvspan(534, 562, zorder=0, alpha=0.2, label='TM 04')
    ax[1].axvspan(572, 598, zorder=0, alpha=0.2, label='TM 05')
    ax[1].axvspan(629, 664, zorder=0, alpha=0.2, label='TM 06')
    ax[1].axvspan(692, 713, zorder=0, alpha=0.2, label='TM 07')
    ax[1].axvspan(714, 738, zorder=0, alpha=0.2, label='TM 08')
    ax[1].axvspan(753, 775, zorder=0, alpha=0.2, label='TM 09')
    ax[1].axvspan(854, 881, zorder=0, alpha=0.2, label='TM 10')

plt.legend()
    # for i,c in zip(np.arange(len(pore_label)),color_pore):
    #     plt.plot((0,0),
    #          (pore_cutoff[2*i],pore_cutoff[2*i+1]),
    #          color=c,
    #          lw=10,
    #          label=pore_label[i],
    #          )

    # ax[1].set_title("SYSTEM %s"%(ax[0]+1))

# ### MISCELLANEOUS ###
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,7.5)
plt.locator_params(axis='y', nbins=6)
plt.tight_layout()
# plt.savefig("%s"%(ifile[:-3]))
plt.show()
