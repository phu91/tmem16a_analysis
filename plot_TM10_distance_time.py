import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
import matplotlib.ticker as ticker  

sns.set_context('notebook')
###############
# FUNCTIONS


def skipping(DATA):
    print("===> SKIP: %s\n"%(skip))
    if skip>1:
        frames = len(DATA.drop_duplicates(['FRAME']))
        skip_period = np.arange(0,frames,skip,dtype='int')
        DATA = DATA[DATA['FRAME'].isin(skip_period)]
        DATA.to_csv("SHORTED_SKIPPED_%s_%s"%(skip,ifile))
        nameout = "SHORTED_SKIPPED_%s_%s"%(skip,ifile)
        print("===> TOTAL FRAMES: %s\n"%(frames))
        print("===> USED FRAMES: %s\n"%(skip_period))
        print("===> NEW INPUT: %s"%(nameout))

        with open("SHORTED_SKIPPED_%s_%s"%(skip,ifile),'a') as file:
            file.write("#SKIPPED %s\n"%(skip))
            file.write("#CUT-OFF %s\n"%(cutoff))
        DATA=pd.read_csv("SHORTED_SKIPPED_%s_%s"%(skip,ifile))
    elif skip ==1:
        DATA=DATA
    return DATA

###############
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT FILE')
parser.add_argument('--cutoff', type=int, default='5',
                    help='Cutoff value between two Centers of Mass. Default 5 (A)')
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
    # delim_whitespace=True,
    comment='#',
    header=0,
    # index_col=0,
    # names=['FRAME','RESNAME_10','RESID_10','CHAIN_10','RESNAME','RESID','CHAIN','DISTANCE'],
    # chunksize=10000
    )
    data = data[data["FRAME"].str.contains("FRAME") == False]
    # print(data)
    data3=data.rename(columns={'TM10': "RESIDUE"})
    # # print(data3)
    data3_melted=data3.melt(["RESIDUE","FRAME"]).rename(columns={'variable':"TM10",'value':"DISTANCE"})
    data3_melted['DISTANCE'] = pd.to_numeric(data3_melted['DISTANCE'])
    # print(data3_melted)
    # plot_data=data3_melted.loc[data3_melted['DISTANCE']<=cutoff]
    # print(plot_data.RESIDUE.str[3:])
    # u=plot_data.sort_values(by='RESIDUE', key=lambda x: x.str[3:])
    # print(plot_data)
    # print(u.sort_values(['FRAME']))
    data3_melted['COUNT']=0
    data3_melted.loc[data3_melted['DISTANCE']<=cutoff,'COUNT']=1
    selected=data3_melted.loc[data3_melted.COUNT==1]
    selected=selected[['RESIDUE','FRAME','COUNT']]
    print(selected)
    print(selected.groupby(['RESIDUE','FRAME'],sort=False).sum())
    # selected=data3_melted.loc[data3_melted['DISTANCE']<=cutoff].sort_values(['FRAME'])
    # print(selected.value_counts().reset_index().rename(columns={0:"COUNT"}).sort_values(['FRAME']))

#     data['RESIDUE']=data['RESNAME']+data['RESID'].astype('str')
#     data['RESIDUE_10']=data['RESNAME_10']+data['RESID_10'].astype('str')
#     data=data[['FRAME','RESIDUE_10','CHAIN_10','RESIDUE','DISTANCE']]
#     data= data.round({"DISTANCE":2})
#     data.to_csv("SHORTED_%s"%(ifile))
#     del data
#     data = pd.read_csv("SHORTED_%s"%(ifile))
#     data = skipping(data)
# else:
#     print("===> Used OLD data!!\n")
#     data = pd.read_csv("SHORTED_%s"%(ifile))
#     data = skipping(data)


# chaina,vmaxa = pivot_data(chain='A')
# chainb,vmaxb = pivot_data(chain='B')

# fig,axes  = plt.subplots(1,2,sharey=False,sharex=True)
# cmap='plasma'
# data_list =[chaina,chainb]
# chain_list=['A','B']

# if vmaxa>vmaxb:
#     vmax=vmaxa
# else:
#     vmax=vmaxb

# for ind,ax in enumerate(axes):
#     g1 = sns.heatmap(data=data_list[ind],
#     cmap=cmap,
#     vmin=0,
#     vmax=vmax,
#     ax=ax,
#     # xticklabels=True, 
#     # yticklabels=True,
#     cbar_kws=({'label':"Number of Contacts",'ticks':np.arange(0, vmax+1)})
#     )
#     g1.set_title("Chain %s | Cut-off: %s (A)"%(chain_list[ind],cutoff))
#     ax.yaxis.set_tick_params(labelsize=8,)
#     ax.yaxis.set_minor_locator(AutoMinorLocator())

# # # # # # ######################
# # # # # ##### MISCELLANEOUS ###
# # # plt.yticks(fontsize=5)
# # # plt.ylim([0,1.1])
# # plt.title("%s |ss%s (A)"%(ifile[:-3],cutoff),va='top')
# plt.suptitle("%s"%(ifile),va='top')
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5*2,3.5*2)   ## Wide x Height
# # # plt.legend()

# plt.tight_layout()
# plt.savefig("%s"%(ifile))
# # plt.show()

