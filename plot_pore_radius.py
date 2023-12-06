import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse,os
from matplotlib import rc
from statannot import add_stat_annotation

###############
# FUNCTIONS

###############
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT FILE REP 1')

args = parser.parse_args()
ifile =  args.input

######################
data = pd.read_csv(ifile,
delim_whitespace=True,
comment='#',
names=['FRAME','CHAIN','ZPOS','RADII'])

# print(data)

os.system("grep COM %s | awk '{print $3,$4,$5}'> COM_%s"%(ifile,ifile))

com_df = pd.read_csv("COM_%s"%(ifile),
delimiter=' ',
names=['CHAIN','FRAME','ZCOM'])

com_average_df = com_df.groupby('CHAIN').mean()
com_a_average = com_average_df.ZCOM[0]
com_b_average = com_average_df.ZCOM[1]
offset = 3

fig, axes = plt.subplots(1,2,sharey=True)

g1 = sns.lineplot(data=data.query("CHAIN=='A'"),
x='ZPOS',
y='RADII',
hue='FRAME',
marker='o',
legend=False,
ax=axes[0])
g1.set_xlim([com_a_average-3,com_a_average+3])
g1.set_ylim([0,5])
g1.axvline(x=com_a_average,linestyle='--',color='b')
g1.axhline(y=1.81,linestyle='--',color='orange')

g2 = sns.lineplot(data=data.query("CHAIN=='B'"),
x='ZPOS',
y='RADII',
hue='FRAME',
marker='s',
ax=axes[1])
g2.set_xlim([com_b_average-3,com_b_average+3])
g2.set_ylim([0,5])
g2.axvline(x=com_b_average,linestyle='-.',color='b',label='COM of PORE')
g2.axhline(y=1.81,linestyle='--',color='orange',label='Chloride radii')


# # ### MISCELLANEOUS ###
plt.legend()
plt.suptitle("%s"%(ifile[:-3]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(25,3)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
# plt.savefig("KDE%s"%(ifile[:-3]))
plt.show()

