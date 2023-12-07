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

def histogram_radius(bins,chain):
    selected_range = np.linspace(left_bound_a,right_bound_a,num=bins,endpoint=True)
    with open("PORE_PROFILE_HISTOGRAM_CHAIN_%s_%s.dat"%(chain,ifile[:-4]),'w+') as ofile:
        for i in range(bins-1):
            in_range=[]
            # for j in data.loc[data.CHAIN=='%s'%(chain)].iterrows():
            #     # print(j[1][3])
            #     if selected_range[i]<=j[1][2]<selected_range[i+1]:
            radius = data.loc[(data.CHAIN=='%s'%(chain)) & (data.ZPOS<=selected_range[i+1]) & (data.ZPOS>=selected_range[i])]
            print(radius.RADII)
    #         in_range.append(j[1][3])
    #         zpos_mid = (selected_range[i]+selected_range[i+1])/2
    #         ofile.write("%s\t%s\t%s\n"%(zpos_mid,np.mean(in_range),np.std(in_range)))
    #         ofile.flush()
    #         # ofile.write("%s %s %s\n"%((selected_range[i]+selected_range[i+1])/2,np.mean(in_range),np.std(in_range)))
    # plot_data = pd.read_csv("PORE_PROFILE_HISTOGRAM_CHAIN_%s_%s.dat"%(chain,ifile[:-4]),delim_whitespace=True,names=['ZPOS','RADII','STD'])
    # plot_data['CHAIN']=chain
    # os.system("rm %s"%("PORE_PROFILE_HISTOGRAM_CHAIN_%s_%s.dat"%(chain,ifile[:-4])))
    # return plot_data
###############
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT FILE REP 1')

args = parser.parse_args()
ifile =  args.input

######################
data = pd.read_csv(ifile,
delim_whitespace=True,
comment='#',
names=['FRAME','CHAIN','ZPOS','RADII'])

os.system("grep COM %s | awk '{print $3,$4,$5}'> COM_%s"%(ifile,ifile))

com_df = pd.read_csv("COM_%s"%(ifile),
delimiter=' ',
names=['CHAIN','FRAME','ZCOM'])

com_average_df = com_df.groupby('CHAIN').mean()
# print(com_average_df)
com_a_average = com_average_df.ZCOM[0]
com_b_average = com_average_df.ZCOM[1]
offset = 3

right_bound_a = com_a_average+offset
left_bound_a = com_a_average-offset
right_bound_b = com_b_average+offset
left_bound_b = com_b_average-offset

#### HISTOGRAM
nbin=5
chain_a = histogram_radius(bins=nbin,chain='A')
chain_b = histogram_radius(bins=nbin,chain='B')

data_2 = pd.concat([chain_a,chain_b])
# print(data_2)

# plt.scatter(chain_a.ZPOS,chain_a.RADII,'o')
# plt.scatter(chain_b.ZPOS,chain_b.RADII,'s')
plt.errorbar(chain_a.ZPOS,chain_a.RADII,chain_a.STD)
plt.errorbar(chain_b.ZPOS,chain_b.RADII,chain_b.STD)

plt.axvline(x=com_a_average,linestyle='-.',color='blue',label='COM of PORE A')
plt.axvline(x=com_b_average,linestyle='-.',color='orange',label='COM of PORE B')
plt.axhline(y=1.81,linestyle='--',color='black',label='Chloride radii')

# fig, axes = plt.subplots(1,2,sharex=True,sharey=True)

# n_open_a  = len(data.loc[(data.CHAIN=='A') & (data.ZPOS<=right_bound_a) & (data.ZPOS>=left_bound_a) & (data.RADII>1.81)].groupby(['FRAME']).mean())
# n_close_a = len(data.loc[(data.CHAIN=='A') & (data.ZPOS<=right_bound_a) & (data.ZPOS>=left_bound_a) & (data.RADII<1.81)].groupby(['FRAME']).mean())
# print(n_open_a)
# n_open_b  = len(data.loc[(data.CHAIN=='B') & (data.ZPOS<=right_bound_b) & (data.ZPOS>=left_bound_b) & (data.RADII>1.81)].groupby(['FRAME']).mean())
# n_close_b = len(data.loc[(data.CHAIN=='B') & (data.ZPOS<=right_bound_b) & (data.ZPOS>=left_bound_b) & (data.RADII<1.81)].groupby(['FRAME']).mean())
# print(n_open_b)

# sns.histplot(data=data.query("CHAIN=='A' & RADII <5"),
# # x='ZPOS',
# x='RADII',
# bins=50,
# ax=axes[0])
# sns.histplot(data=data.query("CHAIN=='B' & RADII <5"),
# # x='ZPOS',
# x='RADII',
# bins=50,
# ax=axes[1])


# g1 = sns.scatterplot(data=data.query("CHAIN=='A' & RADII <3"),
# x='ZPOS',
# y='RADII',
# hue='FRAME',
# marker='o',
# legend=False,
# estimator=np.mean,
# n_boot=10,
# ax=axes[0])
# g1.set_xlim([com_a_average-3,com_a_average+3])
# g1.set_ylim([0,5])
# g1.axvline(x=com_a_average,linestyle='--',color='b')
# g1.axhline(y=1.81,linestyle='--',color='orange')

# g2 = sns.scatterplot(data=data.query("CHAIN=='B' & RADII <3"),
# x='ZPOS',
# y='RADII',
# hue='FRAME',
# estimator=np.mean,
# marker='s',
# ax=axes[1])
# g2.set_xlim([com_b_average-3,com_b_average+3])
# g2.set_ylim([0,5])
# g2.axvline(x=com_b_average,linestyle='-.',color='b',label='COM of PORE')
# g2.axhline(y=1.81,linestyle='--',color='orange',label='Chloride radii')


# # ### MISCELLANEOUS ###
plt.legend()
plt.suptitle("%s"%(ifile[:-3]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(25,3)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
# plt.savefig("KDE%s"%(ifile[:-3]))
plt.show()

