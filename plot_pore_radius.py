import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse,os
from matplotlib import rc
from statannot import add_stat_annotation
import scipy.stats
###############
# FUNCTIONS

def histogram_radius(bins,chain,left_limit,right_limit):
    selected_range = np.linspace(left_limit,right_limit,num=bins+1)
    # check_file = os.path.isfile("HISTOGRAM_CHAIN_%s_%s.dat"%(chain,ifile[:-4]))
    if RW=='Yes':
        with open("HISTOGRAM_CHAIN_%s_%s.dat"%(chain,ifile[:-4]),'w+') as ofile:
            for i in range(bins):
                # for j in data.loc[data.CHAIN=='%s'%(chain)].iterrows():
                #     # print(j[1][3])
                #     if selected_range[i]<=j[1][2]<selected_range[i+1]:
                radii = data.loc[(data.CHAIN=='%s'%(chain)) & (data.ZPOS<=selected_range[i+1]) & (data.ZPOS>=selected_range[i])]
                group_frame = radii.RADII.values
                group_mean=np.mean(group_frame)
                group_standard_dev = np.std(group_frame)
                group_zpos = (selected_range[i]+selected_range[i+1])/2
                # print(group_zpos,group_mean,group_standard_dev)
                ofile.write("%s\t%s\t%s\n"%(group_zpos,group_mean,group_standard_dev))
                ofile.flush()

    plot_data = pd.read_csv("HISTOGRAM_CHAIN_%s_%s.dat"%(chain,ifile[:-4]),delim_whitespace=True,names=['ZPOS','RADII','STD'])
    plot_data['CHAIN']=chain
    #os.system("rm %s"%("PORE_PROFILE_HISTOGRAM_CHAIN_%s_%s.dat"%(chain,ifile[:-4])))
    # print(plot_data)
    return plot_data
###############
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT FILE REP 1')
parser.add_argument('--rw', type=str, default='NO',
                    help='Rewrite data yes/no?. Default is no')
parser.add_argument('--offset', type=int, default='2',
                    help='Offset from COM of PORE')
parser.add_argument('--bins', type=int, default='50',
                    help='Number of Bins. Default 50')
args = parser.parse_args()
ifile =  args.input
RW = args.rw
offset=args.offset
numberOfbins=args.bins

RW = RW.capitalize()
# print(RW)
######################
data = pd.read_csv(ifile,
delim_whitespace=True,
comment='#',
names=['FRAME','CHAIN','ZPOS','RADII'])

# check_file = os.path.isfile("COM_%s"%(ifile))
# print(check_file)
if RW=='Yes':
    print("\n===> Rewriting Data: %s\n"%(RW))
    os.system("grep COM %s | awk '{print $3,$4,$5}'> COM_%s"%(ifile,ifile))
else:
    print("\n===> Rewriting Data: %s\n"%(RW))


com_df = pd.read_csv("COM_%s"%(ifile),
delimiter=' ',
names=['CHAIN','FRAME','ZCOM'])

com_average_df = com_df.groupby('CHAIN').mean()
# print(com_average_df)
com_a_average = com_average_df.ZCOM[0]
com_b_average = com_average_df.ZCOM[1]
# print(com_a_average)
offset = offset

right_bound_a = com_a_average+offset
left_bound_a = com_a_average-offset
right_bound_b = com_b_average+offset
left_bound_b = com_b_average-offset

#### HISTOGRAM
nbin=numberOfbins
chain_a = histogram_radius(bins=nbin,chain='A',left_limit=left_bound_a,right_limit=right_bound_a)
chain_b = histogram_radius(bins=nbin,chain='B',left_limit=left_bound_b,right_limit=right_bound_b)
# print(chain_a)
fig,axes = plt.subplots(1,2,sharey=True)
axes[0].plot(chain_a.ZPOS,chain_a.RADII,color='Tab:blue',marker='o',label='')
axes[0].fill_between(chain_a.ZPOS, chain_a.RADII-chain_a.STD, chain_a.RADII+chain_a.STD, alpha=.25, linewidth=0,color='Tab:blue')
axes[1].plot(chain_b.ZPOS,chain_b.RADII,color='Tab:orange',marker='s',label='')
axes[1].fill_between(chain_b.ZPOS, chain_b.RADII-chain_b.STD, chain_b.RADII+chain_b.STD, alpha=.25, linewidth=0,color='Tab:orange')

axes[0].axvline(x=com_a_average,linestyle=':',color='black',label='ZCOM of PORE')
axes[0].axhline(y=1.81,linestyle='-.',color='red',label='Chloride radius')
axes[0].axhline(y=3.2,linestyle='--',color='blue',label='Hydrated Chloride radius')

axes[1].axvline(x=com_b_average,linestyle=':',color='black',label='ZCOM of PORE')
axes[1].axhline(y=1.81,linestyle='-.',color='red',label='Chloride radius')
axes[1].axhline(y=3.2,linestyle='--',color='blue',label='Hydrated Chloride radius')


axes[0].set_xlim([com_a_average-1,com_a_average+1])
axes[1].set_xlim([com_b_average-1,com_b_average+1])

axes[0].set_title("Chain A")
axes[1].set_title("Chain B")

axes[0].set_xlabel("Z Coordinate")
axes[0].set_ylabel("Pore Radius (A)")
axes[0].set_ylim([0,5])

axes[1].set_xlabel("Z Coordinate")
axes[1].set_ylabel("Pore Radius (A)")
axes[1].set_ylim([0,5])
# # # ### MISCELLANEOUS ###
plt.legend()
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(25,3)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
plt.savefig("%s.png"%(ifile[:-4]))
plt.show()

