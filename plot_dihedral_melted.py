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
                   names=['FRAME','CHAIN','RESIDUE','RESID','DIHEDRAL','ANGLE','SYSTEM','REP'])

res_list=['879','884']
dihedral_list=['CHI1','CHI2']
chain_list=['A','B']

#############################
fig,axes = plt.subplots(2,2,sharex=True,sharey=False,layout="constrained")
# axes = axes.flatten()
data = data[(data["DIHEDRAL"] =='CHI1') | (data["DIHEDRAL"]=='CHI2')]
print(data)

for row in range(len(chain_list)):
    for col in range(len(res_list)):
        g = sns.histplot(
            data=data.query("RESID==%s and CHAIN=='%s'"%(res_list[col],chain_list[row])),
            x="ANGLE",
            hue="DIHEDRAL",
            multiple="stack",
            palette="light:m_r",
            edgecolor=".3",
            linewidth=.5,
            ax=axes[row][col]
            # log_scale=True,
        )
        g.set_title("RES %s | CHAIN %s"%(res_list[col],chain_list[row]))

# ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,6)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
# plt.tight_layout()
plt.savefig("STACKED_%s"%(ifile[:-3]))
plt.show()
