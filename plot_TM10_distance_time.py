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
    delim_whitespace=True,
    comment='#',
    header=0,
    )

    data['NUMBER OF CONTACTS']=1
    u = data.groupby(['FRAME','TM10_CHAIN','OTHER_PART']).sum().reset_index()
    # print(u.sum())
    fig,axes = plt.subplots(2,1,sharex=True,sharey=True)
    chain_list=['A','B']
    for ind,ax in enumerate(axes):
        g = sns.lineplot(data=u.query("TM10_CHAIN=='%s'"%(chain_list[ind])),
        x='FRAME',
        y='NUMBER OF CONTACTS',
        # style='TM10_CHAIN',
        hue='OTHER_PART',
        # markers=True,
        ax=ax
        )
        g.set_title("TM10 CHAIN %s"%(chain_list[ind]))

# # # # # # ######################
# # # # # ##### MISCELLANEOUS ###
# # plt.yticks(fontsize=5)
# # plt.ylim([0,1.1])
# plt.title("%s |ss%s (A)"%(ifile[:-3],cutoff),va='top')
plt.suptitle("%s"%(ifile),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
SCALE=2
plt.gcf().set_size_inches(7.5*SCALE,5.5*SCALE)   ## Wide x Height
# # # plt.legend()
plt.tight_layout()
# plt.savefig("%s"%(ifile))
plt.show()

