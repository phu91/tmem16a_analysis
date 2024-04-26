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

parser.add_argument('--tm46', type=str, default='',
                    help='INPUT Chloride profile')
parser.add_argument('--tm35', type=str, default='',
                    help='INPUT Chloride profile')

parser.add_argument('--factor', type=float, default='1',
                    help='Conversion factor from FRAME to NS. Default 1 frame/ns.'
                    )

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='INPUT PORE Topology')

args = parser.parse_args()

ifile1 =  args.tm46
ifile2 =  args.tm35
systemName = args.system

data1= pd.read_csv(ifile1,comment='#',
delim_whitespace=True,header=0,
)

data2= pd.read_csv(ifile2,comment='#',
delim_whitespace=True,header=0,
)

# print(data1)
# print(data2)

tm46_A = data1.query("chain=='A'").distance.reset_index(drop=True)
tm35_A = data2.query("chain=='A' and resid==567").distance.reset_index(drop=True)
tm46_B = data1.query("chain=='B'").distance.reset_index(drop=True)
tm35_B = data2.query("chain=='B' and resid==567").distance.reset_index(drop=True)

# print(tm46_A)
chainA_dict = {"tm46":tm46_A,
            "tm35":tm35_A}
chainB_dict = {"tm46":tm46_B,
            "tm35":tm35_B}

chainA_df = pd.DataFrame(chainA_dict)
chainB_df = pd.DataFrame(chainB_dict)

# print(chainA_df)

g1=sns.jointplot(data=chainA_df,
            x='tm35',
            y='tm46',
            kind='reg',
            )

g2=sns.jointplot(data=chainB_df,
            x='tm35',
            y='tm46',
            kind='reg',
            )

g1.fig.suptitle("CHAIN A %s"%(systemName))
g2.fig.suptitle("CHAIN B %s"%(systemName))

# data['time (ns)']=data['frame']/5
# grid = plt.GridSpec(3, 2)  #, wspace=0.4, hspace=0.3)
# time_plot_A=plt.subplot(grid[0, 0:])
# time_plot_B=plt.subplot(grid[1, 0:])
# hist_A_plot=plt.subplot(grid[2, 0])
# hist_B_plot=plt.subplot(grid[2, 1])

# # grid.update(hspace=0)

# sns.lineplot(data=data.query("resid==451"),
#             x='time (ns)',
#             y='distance',
#             hue='chain',
#             # style='chain',
#             # markers=True,
#             dashes=False,
#             palette='tab10',
#             ax=time_plot_A)

# sns.lineplot(data=data.query("resid==567"),
#             x='time (ns)',
#             y='distance',
#             hue='chain',
#             # style='chain',
#             # markers=True,
#             dashes=False,
#             palette='tab10',
#             legend=False,
#             ax=time_plot_B)
            
# sns.histplot(data=data.query("resid==451"),
#             x='distance',
#             hue='chain',
#             stat='frequency',
#             bins=20,
#             binrange=[0,30],
#             # y='distance',
#             # color='#1f77b4',
#             multiple='dodge',
#             ax=hist_A_plot)

# sns.histplot(data=data.query("resid==567"),
#             x='distance',
#             hue='chain',
#             stat='frequency',
#             bins=20,
#             binrange=[0,30],
#             # y='distance',
#             multiple='dodge',
#             # color='#1f77b4',
#             ax=hist_B_plot)

# time_plot_A.set_ylim([0,30])
# time_plot_B.set_ylim([0,30])

# # hist_A_plot.set_ylim([0,1])
# # hist_B_plot.set_ylim([0,1])

# hist_A_plot.set_xlim([0,30])
# hist_B_plot.set_xlim([0,30])

# weight='bold'
# time_plot_A.set_title("R451 - K887",fontweight=weight)
# time_plot_B.set_title("K567 - K887",fontweight=weight)
# hist_A_plot.set_title("R451 - K887",fontweight=weight)
# hist_B_plot.set_title("K567 - K887",fontweight=weight)
# # hist_A_plot.set_xlim([0,20])
# # hist_B_plot.set_xlim([0,20])

# time_plot_A.xaxis.set_major_locator(ticker.AutoLocator())
# time_plot_B.xaxis.set_major_locator(ticker.AutoLocator())

# plt.suptitle("%s"%(ifile1))
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,7.5)
# time_plot_A.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
# plt.tight_layout()
g1.fig.savefig("%s_A_CORR.png"%(systemName),dpi=700)
g2.fig.savefig("%s_B_CORR.png"%(systemName),dpi=700)
# plt.show()