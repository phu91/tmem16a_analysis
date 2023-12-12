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
parser.add_argument('--bcolor', type=str, default='bright',
                    help='Background color')
parser.add_argument('--rw', type=str, default='no',
                    help='Rewrite Data. Default: no')
args = parser.parse_args()

ifile =  args.input
background_color = args.bcolor
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
    # chunksize=1000
    )
    # u = data.loc[(data.CHAIN_10=='A') & (data.RESNAME_10=='ASP') & (data.RESID_10==916) & (data.DISTANCE<=cutoff)]
    # print(u.reset_index())
    # Define HELIX UPDATED BASED ON BENDIX
    tm_01 = [327,360]
    tm_02 = [399,439]
    tm_03 = [478,518]
    tm_04 = [534,562]
    tm_05 = [572,598]
    tm_06 = [629,664]
    tm_07 = [692,713]
    tm_08 = [718,740]
    tm_09 = [753,780]
    tm_10 = [887,927]
    
    tm_color_list=sns.color_palette(n_colors=10)
    # print(tm_color_list)
    tm_range_list=[tm_01,tm_02,tm_03,tm_04,tm_05,tm_06,tm_07,tm_08,tm_09,tm_10]
    
    data['TM']='0'
    
    data.loc[(data.RESID<=tm_01[1])&(data.RESID>=tm_01[0]),'TM']='1'
    data.loc[(data.RESID<=tm_02[1])&(data.RESID>=tm_02[0]),'TM']='2'
    data.loc[(data.RESID<=tm_03[1])&(data.RESID>=tm_03[0]),'TM']='3'
    data.loc[(data.RESID<=tm_04[1])&(data.RESID>=tm_04[0]),'TM']='4'
    data.loc[(data.RESID<=tm_05[1])&(data.RESID>=tm_05[0]),'TM']='5'
    data.loc[(data.RESID<=tm_06[1])&(data.RESID>=tm_06[0]),'TM']='6'
    data.loc[(data.RESID<=tm_07[1])&(data.RESID>=tm_07[0]),'TM']='7'
    data.loc[(data.RESID<=tm_08[1])&(data.RESID>=tm_08[0]),'TM']='8'
    data.loc[(data.RESID<=tm_09[1])&(data.RESID>=tm_09[0]),'TM']='9'
    data.loc[(data.RESID<=tm_10[1])&(data.RESID>=tm_10[0]),'TM']='10'
    
    # print(data.loc[data.TM=='2'])
    
    data_sum_a = test(cutoff,'A').sort_values(by="RESID_10",)
    data_sum_b = test(cutoff,'B').sort_values(by="RESID_10",)
    
    # print(data_sum_a)
    data_sum_a.to_csv("TMP_A")
    data_sum_b.to_csv("TMP_B")
    chain_a = pd.read_csv("TMP_A",header=0)
    chain_b = pd.read_csv("TMP_B",header=0)
else:
    print("\n===> Data are NOT overwritten!!\n")
    chain_a = pd.read_csv("TMP_A",header=0)
    chain_b = pd.read_csv("TMP_B",header=0)
    # print(chain_a)

chain_a=chain_a.loc[chain_a.TM!=0]
chain_b=chain_b.loc[chain_b.TM!=0]


chain_a=chain_a.round({"DISTANCE":1})
chain_b=chain_b.round({"DISTANCE":1})
# print(chain_a)

if len(chain_a) == 0 and len(chain_b) == 0:
    print("NO CONTACT FOUND\n")
else:
    markersize=200
    alpha=0.70
    palette='Spectral'

    fig,axes  = plt.subplots(2,1,sharey=True,sharex=True)
    chain_a['CONTACTS']=chain_a['CONTACTS']/(41*855)  ## Maximum number of contacts
    chain_b['CONTACTS']=chain_b['CONTACTS']/(41*855)  ## Maximum number of contacts
    sns.scatterplot(data=chain_a,
    x='RESID_10',
    y='CONTACTS',
    hue='TM',
    hue_order=['1','2','3','4','5','6','7','8','9','10'],
    # style='TM',
    markers=True,
    ax=axes[0],
    # s=markersize,
    alpha=alpha,
    palette=palette,
    size='DISTANCE',
    sizes=(40,200),
    )
    axes[0].set_ylabel("FRACTION OF CONTACTS")
    axes[0].set_xlabel("RESID TM 10 CHAIN A")

    sns.scatterplot(data=chain_b,
    x='RESID_10',
    y='CONTACTS',
    hue='TM',
    hue_order=['1','2','3','4','5','6','7','8','9','10'],
    # style='DISTANCE',
    markers=True,
    ax=axes[1],
    # s=markersize,
    alpha=alpha,
    palette=palette,
    size='DISTANCE',
    sizes=(40,200),
    )
    axes[1].set_ylabel("FRACTION OF CONTACTS")
    axes[1].set_xlabel("RESID TM 10 CHAIN B")

######################

# tm10a_resname, other_resname_a, tm10a_mean_pivot = pivot_data(cutoff=cutoff,chain='A')
# tm10b_resname, other_resname_b, tm10b_mean_pivot = pivot_data(cutoff=cutoff,chain='B')

# # ######################
# fig,axes  = plt.subplots(2,1,sharex=True)
# data_list = [tm10a_mean_pivot,tm10b_mean_pivot]
# chain_list_10= ['A','B']
# chain_list= ['B','A']
# other_resname_list =[other_resname_a,other_resname_b]

# axes =axes.flatten()
# for ind,ax in enumerate(axes):
#     # print(data_list[ind])
#     sns.heatmap(data=data_list[ind],
#     ax=ax,
#     cbar_kws={'label': 'Distance'},
#     cmap='RdPu',
#     vmin=0,
#     vmax=cutoff,
#     xticklabels=True, yticklabels=True,
#     # square=True,
#     # linewidth=.5
#     )
#     ax.set_xlabel("CHAIN %s"%(chain_list[ind]))
#     ax.set_ylabel("TM 10 CHAIN %s"%(chain_list_10[ind]))
#     # v = ax.get_xticklabels()
#     # print(v)
# # ######################
# # # import plotly.express as px
# # # fig = px.imshow(tm10a_mean_pivot,
# #                 # labels=dict(x="Day of Week", y="Time of Day", color="Productivity"),
# #                 # x=['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'],
# #                 # y=['Morning', 'Afternoon', 'Evening']
# #             #    )
# # # fig.update_xaxes(side="top")
# # # fig.show()

# # ######################
# # ##### MISCELLANEOUS ###
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
plt.legend()
plt.tight_layout()
plt.savefig("%s"%(ifile))
#plt.show()

