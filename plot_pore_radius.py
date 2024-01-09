import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse,os
from matplotlib import rc

###############
# FUNCTIONS
sns.set_context('notebook')
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
                    help='Extension factor in Z direction of the PORE')
parser.add_argument('--bins', type=int, default='20',
                    help='Number of Bins. Default 20')
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
com_std_df = com_df.groupby('CHAIN').std()
# print(com_std_df)
com_a_average = com_average_df.ZCOM[0]
com_b_average = com_average_df.ZCOM[1]
com_a_std = com_std_df.ZCOM[0]
com_b_std = com_std_df.ZCOM[1]

nbin=numberOfbins

#### HISTOGRAM
left_bound_a_small = com_a_average-com_a_std*offset
right_bound_a_small = 5
left_bound_b_small = com_b_average-com_b_std*offset
right_bound_b_small = 5

chain_a_small = histogram_radius(bins=nbin,chain='A',left_limit=left_bound_a_small,right_limit=right_bound_a_small)
chain_b_small = histogram_radius(bins=nbin,chain='B',left_limit=left_bound_b_small,right_limit=right_bound_b_small)
# print(chain_a.ZPOS)

# Initialize the grid
grid = plt.GridSpec(1, 2, wspace=0.3, hspace=0.3)
# make subplots
g1 = plt.subplot(grid[0])
g2 = plt.subplot(grid[1])

color_a='violet'
color_b='deepskyblue'
alpha=0.25
g1.axvline(x=com_a_average,linestyle='--',lw=3,color='black',label='Center of Mass of Pore A')
g1.axhline(y=1.81,linestyle='-.',lw=3,color='coral',label='Hydrated Cl-')
g1.axhline(y=3.2,linestyle='-.',lw=3,color='blue',label='Cl-')

g2.axvline(x=com_b_average,linestyle='--',lw=3,color='black',label='Center of Mass of Pore B')
g2.axhline(y=1.81,linestyle='-.',lw=3,color='coral',label='Hydrated Cl-')
g2.axhline(y=3.2,linestyle='-.',lw=3,color='blue',label='Cl-')


g1.plot(chain_a_small.ZPOS,chain_a_small.RADII,color=color_a,marker=None,label='',lw=6,markersize=10)
g1.fill_between(chain_a_small.ZPOS, chain_a_small.RADII-chain_a_small.STD, chain_a_small.RADII+chain_a_small.STD, alpha=alpha, linewidth=0,color=color_a)
g2.plot(chain_b_small.ZPOS,chain_b_small.RADII,color=color_b,marker=None,label='',lw=6,markersize=10)
g2.fill_between(chain_b_small.ZPOS, chain_b_small.RADII-chain_b_small.STD, chain_b_small.RADII+chain_b_small.STD, alpha=alpha, linewidth=0,color=color_b)

g1.set_xlim([com_a_average-com_a_std*offset,5])
g2.set_xlim([com_b_average-com_b_std*offset,5])

fontsize=20
g1.set_title("Chain A",weight='bold',fontsize=fontsize)
g2.set_title("Chain B",weight='bold',fontsize=fontsize)

g1.set_xlabel("Z Coordinate",fontsize=fontsize,weight='bold')
g1.set_ylabel("Radius (Å)",fontsize=fontsize,weight='bold')
g1.set_ylim([0,4])

g2.set_xlabel("Z Coordinate",fontsize=fontsize,weight='bold')
g2.set_ylabel("Radius (Å)",fontsize=fontsize,weight='bold')
# g2.set_xlim([-10,10])
g2.set_ylim([0,4])


# # # ### MISCELLANEOUS ###
g1.legend(loc='lower right')
g2.legend(loc='lower right')

plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
SCALE=1
plt.gcf().set_size_inches(7.5*SCALE,6.5*SCALE)   ## Wide x Height
plt.tight_layout()
plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
# plt.show()

