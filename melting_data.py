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
parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='NAME OF THE SYSTEM')
args = parser.parse_args()

ifile =  args.input
systemname = args.system

data = pd.read_csv(ifile,comment='#',
                   delim_whitespace=True,
                   names=['FRAME','CHAIN','RESIDUE NAME','RESID','OMEGA','PHI','PSI','CHI1','CHI2','SYSTEM'])
                   
data_melt = data.melt(id_vars=['FRAME','CHAIN','RESIDUE NAME','RESID'])
data_melt = data_melt.query("variable!='SYSTEM'")
data_melt['SYSTEM']=systemname
data_melt.to_csv("NEW_%s"%(ifile),index=False,sep='\t',header=['#FRAME','CHAIN','RESIDUE','RESID','DIHEDRAL','ANGLE','SYSTEM'])
# g = sns.kdeplot(data=data,
#             x='CHI1',
#             y='CHI2',
#             style='CHAIN',
#             hue='RESID',
#             # markers=True,
#             palette='Set1',
#             )
# g = sns.scatterplot(data=data,
#             x='CHI1',
#             y='CHI2',
#             style='CHAIN',
#             hue='RESID',
#             # markers=True,
#             palette='Set1',
#             )
# color_list = sns.color_palette('Set2')


# ### MISCELLANEOUS ###
# plt.suptitle("%s"%(ifile[:-4]),va='top')
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,6)   ## Wide x Height
# # plt.locator_params(axis='both', nbins=5)
# plt.tight_layout()
# plt.savefig("%s"%(ifile[:-3]))
# plt.show()
