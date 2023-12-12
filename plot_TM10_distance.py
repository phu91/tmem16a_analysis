import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from statannot import add_stat_annotation

sns.set_context('notebook')
###############
# FUNCTIONS

def pivot_data(cutoff,chain):
    # print(data)
    cutoff_data = data.loc[data.DISTANCE<=cutoff]
    # cutoff_data = data
    u = cutoff_data.groupby(['FRAME','RESNAME_10','RESID_10','CHAIN_10'])
    tm10_define = np.arange(887,928)
    other_chain= np.arange(73,928)
    # print(len(tm10_define),len(other_chain))

    shorten_data_10 = cutoff_data.loc[cutoff_data.CHAIN_10=='%s'%(chain)]
    u = shorten_data_10.groupby(['RESNAME_10','RESID_10','RESNAME','RESID','CHAIN'])
    u_mean = u.mean()
    del u_mean['FRAME']
    u_mean = u_mean.reset_index()
    # print(u_mean)
    resname_10_list = u_mean.RESNAME_10
    resname_list = u_mean.RESNAME
    u_mean_pivot = u_mean.pivot_table(index='RESID_10',columns='RESID',values='DISTANCE',fill_value=cutoff+100)
    # print(u_mean_pivot)
    del u, shorten_data_10
    return resname_10_list,resname_list,u_mean_pivot

###############
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT FILE REP 5')
parser.add_argument('--cutoff', type=int, default='5',
                    help='Cutoff value to define BENDING. Default: 5 A')
parser.add_argument('--bcolor', type=str, default='bright',
                    help='Background color')
args = parser.parse_args()

ifile =  args.input
background_color = args.bcolor
cutoff = args.cutoff

######################
data = pd.read_csv(ifile,
delim_whitespace=True,
comment='#',
names=['FRAME','RESNAME_10','RESID_10','CHAIN_10','RESNAME','RESID','CHAIN','DISTANCE'],
# chunksize=1000
)
######################

tm10a_resname, other_resname_a, tm10a_mean_pivot = pivot_data(cutoff=cutoff,chain='A')
tm10b_resname, other_resname_b, tm10b_mean_pivot = pivot_data(cutoff=cutoff,chain='B')

# ######################
fig,axes  = plt.subplots(2,1,sharex=True)
data_list = [tm10a_mean_pivot,tm10b_mean_pivot]
chain_list_10= ['A','B']
chain_list= ['B','A']
other_resname_list =[other_resname_a,other_resname_b]

axes =axes.flatten()
for ind,ax in enumerate(axes):
    # print(data_list[ind])
    sns.heatmap(data=data_list[ind],
    ax=ax,
    cbar_kws={'label': 'Distance'},
    cmap='RdPu',
    vmin=0,
    vmax=cutoff,
    xticklabels=True, yticklabels=True,
    # square=True,
    # linewidth=.5
    )
    ax.set_xlabel("CHAIN %s"%(chain_list[ind]))
    ax.set_ylabel("TM 10 CHAIN %s"%(chain_list_10[ind]))
    # v = ax.get_xticklabels()
    # print(v)
# ######################
# # import plotly.express as px
# # fig = px.imshow(tm10a_mean_pivot,
#                 # labels=dict(x="Day of Week", y="Time of Day", color="Productivity"),
#                 # x=['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'],
#                 # y=['Morning', 'Afternoon', 'Evening']
#             #    )
# # fig.update_xaxes(side="top")
# # fig.show()

# ######################
# ##### MISCELLANEOUS ###
# plt.xticks(rotation=72)
# # plt.ylim([0,1.1])
# plt.title("%s |ss%s (A)"%(ifile[:-3],cutoff),va='top')
plt.suptitle("%s"%(ifile),va='top')
axes[0].set_title("Cut-off: %s (A)"%(cutoff))
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,10)   ## Wide x Height
# plt.locator_params(axis='y', nbins=89)

plt.tight_layout()
# plt.savefig("%s"%(ifile))
plt.show()

