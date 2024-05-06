import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
import matplotlib.ticker as ticker  
import colorcet as cc

# sns.set_context('talk')
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
data = pd.read_csv(ifile,
delim_whitespace=True,
comment='#',
header=None,
names=['FRAME','TM10_CHAIN', 'TM10','OTHER','OTHER_PART','DISTANCE']
)
data['NUMBER OF CONTACTS']=1
# print(data)
u = data.groupby(['FRAME','TM10_CHAIN','OTHER_PART']).sum().reset_index()
u['Time (ns)']=u['FRAME']*2/10
chain_list=['A','B']
palette = sns.color_palette(cc.glasbey_hv, n_colors=20)
u=u[['TM10_CHAIN','OTHER_PART','NUMBER OF CONTACTS','Time (ns)']]
if len(u.query("OTHER_PART=='TM10'")) ==0:
#     # print("No")
    add_tm10_dummy ={'TM10_CHAIN':['A','B'],
                     'OTHER_PART':['TM10','TM10'],
                     'NUMBER OF CONTACTS':[0,0],
                     'Time (ns)':[0,0]
                     }
    u=pd.concat([u,pd.DataFrame(add_tm10_dummy)])

u = u.sort_values('OTHER_PART')
# g=sns.violinplot(data=u.query("OTHER_PART=='L23' or OTHER_PART=='L45' or OTHER_PART=='Nterm' or OTHER_PART=='TM10' or OTHER_PART=='TM3'"),   
g=sns.violinplot(data=u.query("OTHER_PART=='L23' or OTHER_PART=='Nterm' or OTHER_PART=='TM10'"),   
x='OTHER_PART',
y='NUMBER OF CONTACTS',
hue='TM10_CHAIN',
# hue_order=['Nterm','L23','TM3','TM10'],
split=True,
palette='cool_r'
)
plt.ylim([0,35])
fontsize=25
weight='bold'
plt.tick_params(axis='x', labelsize=15)
plt.xlabel("Other Chain",fontsize=fontsize,weight=weight)
plt.ylabel("Number of Contacts",fontsize=fontsize,weight=weight)


##### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile),va='top',weight='bold')
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.fonttype'] = 42
SCALE=1.
plt.gcf().set_size_inches(7.5*SCALE,5.5*SCALE)   ## Wide x Height
plt.tight_layout()
plt.savefig("%s_KIND.png"%(ifile),dpi=700)
# plt.show()

