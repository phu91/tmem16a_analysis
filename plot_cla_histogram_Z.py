import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
import colorcet as cc
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker

sns.set(style="ticks", context="notebook")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT Chloride profile')

# parser.add_argument('--topol', type=str, default='',
#                     help='INPUT PORE Topology')

parser.add_argument('--factor', type=float, default='1',
                    help='Conversion factor from FRAME to NS. Default 1 frame/ns.'
                    )

args = parser.parse_args()

ifile =  args.input
# pore_top = args.topol
factor = args.factor

data = pd.read_csv(ifile,comment='#',
delim_whitespace=True,
names=['FRAME','INDEX','CHAIN','ZPOS','POSITION_PORE','DCOM_PORE'])
binrange=[-25,25]
data['ns']=data.FRAME*factor
# print(data)
plot = sns.histplot(data=data,
y='ZPOS',
hue='CHAIN',
stat='density',
binwidth=1,
multiple='stack',
# kde=True,
binrange=binrange,
palette='Set2',
common_norm=True,
)

plot.set_ylim(binrange)

# # # # ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,5.5)
# plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()

# plt.savefig("%s.png"%(ifile[:-4]),dpi=700)

plt.show()