import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
import matplotlib.ticker as ticker  
import colorcet as cc

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
    header=None,
    names=['FRAME','TM10_CHAIN', 'TM10','OTHER','OTHER_PART','DISTANCE']
    )

    data['NUMBER OF CONTACTS']=1
    # print(data)
    u = data.groupby(['FRAME','TM10_CHAIN','OTHER_PART']).sum().reset_index()
    # print(u)
    u['Time (ns)']=u['FRAME']*2/10
    fig,axes = plt.subplots(2,1,sharex=True,sharey=True)
    chain_list=['A','B']
    palette = sns.color_palette(cc.glasbey_hv, n_colors=20)
    
    for ind,ax in enumerate(axes):
        g = sns.lineplot(data=u.query("TM10_CHAIN=='%s'"%(chain_list[ind])),
        x='Time (ns)',
        y='NUMBER OF CONTACTS',
        # style='TM10_CHAIN',
        hue='OTHER_PART',
        hue_order=['TM1','TM2','TM3','TM4','TM5','TM6','TM7','TM8','TM9','TM10','L12','L23','L34','L45','L56','L67','L78','L89','L910','Nterm'],
        # markers=True,
        palette=palette,
        lw=3,
        alpha=0.95,
        ax=ax
        )
        ax.set_ylim([0,30])
        weight='bold'
        fontsize=30
        g.set_title("TM10 CHAIN %s"%(chain_list[ind]),weight=weight,fontsize=fontsize)
        ax.tick_params(axis='both', labelsize='large',)
        ax.set_xlabel("Time (ns)",weight=weight,fontsize=fontsize)
        ax.set_ylabel("Number of Contacts", weight=weight,fontsize=fontsize)
    axes[1].get_legend().remove()
    legend=axes[0].legend(loc="upper right",
    # bbox_to_anchor=(.5, -0.4), 
    ncol=5, 
    title="Other Chain",
    frameon=False,
    )
    legend.get_title().set_fontweight('bold')
    for line in legend.get_lines():
        line.set_linewidth(0)
        line.set_marker('s')
        line.set_markersize(15)
    for text in legend.get_texts():
        text.set_fontweight("medium")

##### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile),va='top',weight='bold')
plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'weight':'bold'})
plt.rcParams['pdf.fonttype'] = 42
SCALE=1.5
plt.gcf().set_size_inches(7.5*SCALE,7.5*SCALE)   ## Wide x Height
plt.tight_layout()
plt.savefig("%s.png"%(ifile),dpi=700)
# plt.show()

