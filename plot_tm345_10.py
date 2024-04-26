import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker

sns.set(style="ticks", context="notebook")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT Chloride profile')

parser.add_argument('--factor', type=float, default='1',
                    help='Conversion factor from FRAME to NS. Default 1 frame/ns.'
                    )
# parser.add_argument('--system', type=str, default='UNKNOWN',
#                     help='Default UNKNOWN.'
#                     )

args = parser.parse_args()

ifile =  args.input
# systemName = args.system

data = pd.read_csv(ifile,comment='#',
delim_whitespace=True,header=0,
)

data = data.query("resid==567")
data['time (ns)']=data['frame']/5

grid = plt.GridSpec(2, 2)  #, wspace=0.4, hspace=0.3)
time_plot=plt.subplot(grid[0, 0:])
hist_A_plot=plt.subplot(grid[1, 0])
hist_B_plot=plt.subplot(grid[1, 1])

sns.lineplot(data=data,
            x='frame',
            y='distance',
            hue='chain',
            palette='tab10',
            ax=time_plot)

sns.histplot(data=data.query("chain=='A'"),
            x='distance',
            # y='distance',
            color='#1f77b4',
            ax=hist_A_plot)

sns.histplot(data=data.query("chain=='B'"),
            x='distance',
            # y='distance',
            color='#ff7f0e',
            ax=hist_B_plot)

time_plot.set_ylim([0,20])
hist_A_plot.set_xlim([0,20])
hist_B_plot.set_xlim([0,20])

plt.suptitle("%s"%(ifile))
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,5.5)
plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()
plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
# plt.show()