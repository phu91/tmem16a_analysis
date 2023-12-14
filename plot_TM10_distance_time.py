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
parser.add_argument('--cutoff', type=int, default='10',
                    help='Cutoff value between two Centers of Mass. Default 10 (A)')
parser.add_argument('--skip', type=int, default='0',
                    help='[0,Read data fully],[1,Read skipped data],[skip>1,READ NORMAL DATA and GENERATE SKIPPED DATA]')
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

    data['RESIDUE']=data['RESNAME']+data['RESID'].astype('str')
    data['RESIDUE_10']=data['RESNAME_10']+data['RESID_10'].astype('str')
    data=data[['FRAME','RESIDUE_10','CHAIN_10','RESIDUE','DISTANCE']]
    data= data.round({"DISTANCE":2})
    data.to_csv("SHORTED_%s"%(ifile))

else:
    print("===> Used OLD data!!\n")
    if skip==0:  # Read data fully
        data = pd.read_csv("SHORTED_%s"%(ifile))
    elif skip==1: # Read skipped data
        data=pd.read_csv("SHORTED_SKIPPED_%s"%(ifile))
    elif skip>1: #Read data fully and generate skipped data
        data = pd.read_csv("SHORTED_%s"%(ifile))
        frames = len(data.drop_duplicates(['FRAME']))
        print("===> TOTAL FRAMES: %s\n"%(frames))
        skip_period = np.arange(0,frames,skip,dtype='int')
        print("===> USED FRAMES: %s\n"%(skip_period))
        data = data[data['FRAME'].isin(skip_period)]
        data.to_csv("SHORTED_SKIPPED_%s"%(ifile))
        with open("SHORTED_SKIPPED_%s"%(ifile),'a') as file:
            file.write("#SKIPPED %s\n"%(skip))
            file.write("#CUT-OFF %s\n"%(cutoff))

u_a = data.loc[(data.CHAIN_10=='B')]
u_a1 = u_a.groupby(['FRAME','RESIDUE']).mean().reset_index()
u_a1['CONTACT']=0
u_a1.loc[(u_a1.DISTANCE<=cutoff),'CONTACT']=1
u_a1.loc[(u_a1.DISTANCE>cutoff),'CONTACT']=0
u_a2 = pd.pivot_table(u_a1,index='RESIDUE',columns='FRAME',values='CONTACT')

u_b = data.loc[(data.CHAIN_10=='B')]
u_b1 = u_b.groupby(['FRAME','RESIDUE']).mean().reset_index()
u_b1['CONTACT']=0
u_b1.loc[(u_b1.DISTANCE<=cutoff),'CONTACT']=1
u_b1.loc[(u_b1.DISTANCE>cutoff),'CONTACT']=0
u_b2 = pd.pivot_table(u_b1,index='RESIDUE',columns='FRAME',values='CONTACT')

fig,axes  = plt.subplots(1,2,sharey=True,sharex=False)
cmap='RdPu'
data_list =[u_a2,u_b2]
chain_list=['A','B']
for ind,ax in enumerate(axes):
    g1 = sns.heatmap(data=data_list[ind],
    cmap=cmap,
    # vmin=0,
    # vmax=1,
    # robust=True,
    ax=ax,
    # xticklabels=True, 
    # yticklabels=True,
    # cbar_kws=({'label':"Distance (A)"})
    cbar=False
    )
    g1.set_title("TM 10 Chain %s | Cut-off: %s (A)"%(chain_list[ind],cutoff))
    ax.yaxis.set_tick_params(labelsize=8)

    # g1.set_yticklabels(g1.get_yticks(), size = 5)

# # # # ######################
# # # ##### MISCELLANEOUS ###
# plt.yticks(fontsize=5)
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

