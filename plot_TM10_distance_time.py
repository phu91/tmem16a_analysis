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
    # print(cutoff_data)
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

def test(cutoff,chain):
    cutoff_data = data.loc[(data.DISTANCE<=cutoff) & (data.CHAIN==chain)]
    cutoff_data['CONTACTS']=1
    del cutoff_data['FRAME']
    # print(cutoff_data)
    u = cutoff_data.groupby(['RESNAME_10','RESID_10','CHAIN_10'])
    tm10_define = np.arange(887,928)
    other_chain= np.arange(73,928)
    # print(len(tm10_define),len(other_chain))

    shorten_data_10 = cutoff_data
    u = shorten_data_10.groupby(['RESNAME_10','RESID_10','TM'])
    u2 = u.sum().reset_index()
    del u2['RESID']
    del u2['DISTANCE']
    u3 = u.mean().reset_index()
    del u3['RESID']
    del u3['CONTACTS']
    # print(u2)
    # print(u3)
    new_data = u2.merge(u3)

    return new_data

    # return u2
###############
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT FILE')
parser.add_argument('--cutoff', type=int, default='5',
                    help='Cutoff value to define BENDING. Default: 5 A')
parser.add_argument('--skip', type=int, default='1',
                    help='Skip Frames')
parser.add_argument('--rw', type=str, default='no',
                    help='Rewrite Data. Default: no')
args = parser.parse_args()

ifile =  args.input
skip = args.skip
cutoff = args.cutoff
RW = args.rw
# print(RW)
######################
if RW=='yes':
    print("\n===> Data are overwritten!!\n")
    data = pd.read_csv(ifile,
    delim_whitespace=True,
    comment='#',
    names=['FRAME','RESNAME_10','RESID_10','CHAIN_10','RESNAME','RESID','CHAIN','DISTANCE'],
    # chunksize=10000
    )
    frames = len(data.drop_duplicates(['FRAME']))
    print("===> TOTAL FRAMES: %s\n"%(frames))
    skip=skip
    skip_period = np.arange(0,frames,skip,dtype='int')
    print("===> USED FRAMES: %s\n"%(skip_period))
    data = data[data['FRAME'].isin(skip_period)]
    # print("===> NEW DATAFRAME: %s\n"%(frames))
    
    data['RESIDUE']=data['RESNAME']+data['RESID'].astype('str')
    data['RESIDUE_10']=data['RESNAME_10']+data['RESID_10'].astype('str')
    data[['FRAME','RESIDUE_10','RESIDUE','DISTANCE']]
    data.to_csv("SHORTED_%s"%(ifile))
else:
    data = pd.read_csv("SHORTED_%s"%(ifile))

u_a = data.loc[(data.CHAIN_10=='A')]
u_a.loc[(data.DISTANCE>=cutoff),'DISTANCE']=cutoff
u_a1 = u_a.groupby(['FRAME','RESIDUE']).mean().reset_index()
# print(u_a1.loc[u_a1.FRAME==1])
u_a2 = pd.pivot_table(u_a1,index='RESIDUE',columns='FRAME',values='DISTANCE')
u_b = data.loc[(data.CHAIN_10=='B')]
u_b.loc[(data.DISTANCE>=cutoff),'DISTANCE']=cutoff
u_b1 = u_b.groupby(['FRAME','RESIDUE']).mean().reset_index()
u_b2 = pd.pivot_table(u_b1,index='RESIDUE',columns='FRAME',values='DISTANCE')
fig,axes  = plt.subplots(2,1,sharey=True,sharex=True)
cmap='magma_r'
g1 = sns.heatmap(data=u_a2,
cmap=cmap,
vmin=0,
vmax=cutoff,
ax=axes[0],
# xticklabels=True, 
# yticklabels=True,
cbar_kws=({'label':"Distance (A)"})
)
g1.set_title("TM 10 Chain A | Cut-off: %s (A)"%(cutoff))
g2 = sns.heatmap(data=u_b2,
cmap=cmap,
vmin=0,
vmax=cutoff,
ax=axes[1],
# xticklabels=True, 
# yticklabels=True,
cbar_kws=({'label':"Distance (A)"})
)
g2.set_title("TM 10 Chain B | Cut-off: %s (A)"%(cutoff))

# # # # ######################
# # # ##### MISCELLANEOUS ###
plt.xticks(rotation=72)
# plt.ylim([0,1.1])
# plt.title("%s |ss%s (A)"%(ifile[:-3],cutoff),va='top')
plt.suptitle("%s"%(ifile),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,10)   ## Wide x Height
# plt.locator_params(axis='y', nbins=89)
# plt.legend()
plt.tight_layout()
plt.savefig("%s"%(ifile))
# plt.show()

