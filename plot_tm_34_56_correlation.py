import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
import colorcet as cc
from sklearn.metrics import r2_score

sns.set(style="ticks", context="notebook")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT Chloride profile')


args = parser.parse_args()

ifile =  args.input

data = pd.read_csv(ifile,comment='#',
delim_whitespace=True,names=['frame','resid','chain','distance']
)

data['time']=data['frame']*0.2

resid567_A = data.query("chain=='A' and resid==567").distance
resid584_A = data.query("chain=='A' and resid==584").distance
resid567_B = data.query("chain=='B' and resid==567").distance
resid584_B = data.query("chain=='B' and resid==584").distance

fig,axes = plt.subplots(1,2)

axes[0].scatter(x=resid567_A,
y=resid584_A,
color='blue')
z = np.polyfit(resid567_A, resid584_A, 1)
p = np.poly1d(z)
axes[0].plot(resid567_A,p(resid567_A),"r--")
axes[0].set_title("CHAIN A")

axes[1].scatter(x=resid567_B,
y=resid584_B,
color='orange')
z = np.polyfit(resid567_B, resid584_B, 1)
p = np.poly1d(z)
axes[1].plot(resid567_B,p(resid567_B),"r--")
axes[1].set_title("CHAIN B")

plt.suptitle("%s"%(ifile))
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,5.5)
plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()
plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
plt.show()