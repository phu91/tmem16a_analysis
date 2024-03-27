import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
import colorcet as cc
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator, AutoMinorLocator

sns.set(style="ticks", context="notebook")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT Chloride profile')

parser.add_argument('--factor', type=float, default='1.0',
                    help='Time convertion. Default 1 ns/frame')

args = parser.parse_args()

ifile =  args.input
time_factor = args.factor

data = pd.read_csv(ifile,comment='#',
delim_whitespace=True,
names=['frame','resname','resid','chain_a','chain_b'])
data['time']=data['frame']*time_factor
query_data=data.query("637<= resid <=648").dropna().replace(['I','H','T','B','E','G','P','S','~'],[1,0,0,0,0,0,0,0,0])

total_residues = query_data['resid'].max()-query_data['resid'].min()+1

pi_helix_propensity = query_data.groupby(['frame']).sum()/total_residues
print(pi_helix_propensity)
# pi_helix_propensity['frame']=pi_helix_propensity['frame']*0.2

g = sns.lineplot(data=pi_helix_propensity,
                x='time',
                y='chain_a',
                label='CHAIN A'
                )
g = sns.lineplot(data=pi_helix_propensity,
                x='time',
                y='chain_b',
                label='CHAIN B'
                )

g.set_xlabel("Time (ns)")
g.set_ylabel("Pi-Helix Propensity of TM6 Hinge")
g.set_ylim([0,1])
# # # ### MISCELLANEOUS ###
# plt.suptitle("%s"%(ifile[:-4]),va='top')
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,5.5)
# plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
# plt.tight_layout()

plt.savefig("%s.png"%(ifile[:-4]),dpi=700)

# plt.show()
