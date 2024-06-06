import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from statannot import add_stat_annotation

sns.set_context("notebook")
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to RMSF profile')

args = parser.parse_args()

ifile =  args.input
#FRAME L_A L_B theta_A theta_B
data = pd.read_csv(ifile,
                 names=['Frame','l_A','l_B','theta_A','theta_B'],
                 comment='#',
                 delim_whitespace=True)

u = data.std()
print(u.theta_A/u.theta_B)

# ### MISCELLANEOUS ###
# plt.suptitle("%s"%(ifile[:-4]),va='top')
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,10)   ## Wide x Height
# # plt.locator_params(axis='both', nbins=5)
# plt.tight_layout()
# plt.savefig("%s"%(ifile[:-3]))
plt.show()
