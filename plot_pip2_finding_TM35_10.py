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

parser.add_argument('--filter', nargs='+', help='<Required> Set flag', required=True)

args = parser.parse_args()

ifile =  args.input
filter= args.filter

# data = pd.read_csv(ifile,comment='#',
# delim_whitespace=True,names=['frame','resid','position','binding']
# )

# u = data.query("resid==412 and position=='567A_887B'")
# v = data.query("resid==430 and position=='567B_887A'")
# # print(v)

# ax = sns.lineplot(data=u,
# x='frame',
# y='binding',
# marker='o',
# )

# ax = sns.lineplot(data=v,
# x='frame',
# y='binding',
# marker='s',
# )

# plt.suptitle("%s"%(ifile))
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,5.5)
# plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
# plt.tight_layout()
# plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
plt.show()