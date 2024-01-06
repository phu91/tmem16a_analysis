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
                    help='Extension in Z direction of the PORE')
parser.add_argument('--bins', type=int, default='10',
                    help='Number of Bins. Default 10')
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

nbin=numberOfbins

#### HISTOGRAM
right_bound_a_small = com_a_average+offset
left_bound_a_small = com_a_average-offset
right_bound_b_small = com_b_average+offset
left_bound_b_small = com_b_average-offset

chain_a_small = histogram_radius(bins=nbin,chain='A',left_limit=left_bound_a_small,right_limit=right_bound_a_small)
chain_b_small = histogram_radius(bins=nbin,chain='B',left_limit=left_bound_b_small,right_limit=right_bound_b_small)
# print(chain_a.ZPOS)

#### HISTOGRAM
offset2=50
nbin2=nbin+20
right_bound_a = com_a_average+offset2
left_bound_a = com_a_average-offset2
right_bound_b = com_b_average+offset2
left_bound_b = com_b_average-offset2

chain_a = histogram_radius(bins=nbin2,chain='A',left_limit=left_bound_a,right_limit=right_bound_a)
chain_b = histogram_radius(bins=nbin2,chain='B',left_limit=left_bound_b,right_limit=right_bound_b)
# print(chain_a.ZPOS)

# Initialize the grid
grid = plt.GridSpec(2, 2, wspace=0.4, hspace=0.3)
# make subplots
g1 = plt.subplot(grid[0, 0])
g2 = plt.subplot(grid[1, 0])
g3 = plt.subplot(grid[:2,1])

# fig,axes = plt.subplots(1,2,sharey=True)
color_a='violet'
color_b='deepskyblue'
alpha=0.25
g1.axvline(x=com_a_average,linestyle=':',color='black',label='Center of Mass of Pore A')
g1.axhline(y=1.81,linestyle='-.',color='coral',label='')
g1.axhline(y=3.2,linestyle='-.',color='blue',label='')

g2.axvline(x=com_b_average,linestyle='--',color='black',label='Center of Mass of Pore B')
g2.axhline(y=1.81,linestyle='-.',color='coral',label='')
g2.axhline(y=3.2,linestyle='-.',color='blue',label='')

g3.axhline(y=com_a_average,linestyle=':',color='black',label='')
g3.axhline(y=com_b_average,linestyle='--',color='black',label='')
g3.axvline(x=1.81,linestyle='-.',color='coral',label='Chloride radius')
g3.axvline(x=3.2,linestyle='-.',color='blue',label='Hydrated Chloride radius')
#################

g1.plot(chain_a_small.ZPOS,chain_a_small.RADII,color=color_a,marker='o',label='',lw=3,markersize=10)
g1.fill_between(chain_a_small.ZPOS, chain_a_small.RADII-chain_a_small.STD, chain_a_small.RADII+chain_a_small.STD, alpha=alpha, linewidth=0,color=color_a)
g2.plot(chain_b_small.ZPOS,chain_b_small.RADII,color=color_b,marker='s',label='',lw=3,markersize=10)
g2.fill_between(chain_b_small.ZPOS, chain_b_small.RADII-chain_b_small.STD, chain_b_small.RADII+chain_b_small.STD, alpha=alpha, linewidth=0,color=color_b)

g3.fill_betweenx(chain_b.ZPOS, chain_b.RADII+chain_b.STD, chain_b.RADII-chain_b.STD, alpha=alpha, linewidth=0,color=color_b)
g3.fill_betweenx(chain_a.ZPOS, chain_a.RADII+chain_a.STD, chain_a.RADII-chain_a.STD, alpha=alpha, linewidth=0,color=color_a)
g3.plot(chain_b.RADII,chain_b.ZPOS,color=color_b,marker='s',label='',lw=3,markersize=0)
g3.plot(chain_a.RADII,chain_a.ZPOS,color=color_a,marker='o',label='',lw=3,markersize=0)

g1.set_xlim([com_a_average-offset,com_a_average+offset])
g2.set_xlim([com_b_average-offset,com_b_average+offset])

g1.set_title("Chain A",weight='bold')
g2.set_title("Chain B",weight='bold')

g1.set_xlabel("Z Coordinate",weight='bold')
g1.set_ylabel("Radius (Å)",weight='bold')
g1.set_ylim([1,4])

g2.set_xlabel("Z Coordinate",weight='bold')
g2.set_ylabel("Radius (Å)",weight='bold')
g2.set_ylim([1,4])

g3.set_ylabel("Z Coordinate",weight='bold')
g3.set_xlabel("Radius (Å)",weight='bold')
# g3.set_ylim([0,4])

# # # ### MISCELLANEOUS ###
g1.legend()
g2.legend()
g3.legend()
plt.suptitle("%s"%(ifile[:-4]),va='top',weight='bold')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
SCALE=1
plt.gcf().set_size_inches(7.5*SCALE,6.5*SCALE)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
# plt.tight_layout()
plt.savefig("%s.png"%(ifile[:-4]))
# plt.show()

