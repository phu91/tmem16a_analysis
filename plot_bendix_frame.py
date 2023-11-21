import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from statannot import add_stat_annotation

###############
# FUNCTIONS

def find_headers():
    nres=endRes-startRes+1
    headers=['FRAME']
    for i in np.arange(0,nres*2):  # FOR two chain A and B
        if i <nres:
            i = str(i+startRes)+'A'
            headers.append(i)
        else:
            i = str(i+startRes-nres)+'B'
            headers.append(i)
    return headers

def read_data(inputfile):
    df = pd.read_csv(inputfile,delim_whitespace=True,skiprows=13,
                    header=None,names=find_headers())
    df = df.drop_duplicates()
    return df

def fraction_of_angles(dataframe,cutoff):
    nframe = len(dataframe)
    header_noFRAME = find_headers()
    header_noFRAME.remove('FRAME')
    fraction_angles = []
    for i in header_noFRAME:
        u = dataframe.loc[dataframe[i]>cutoff]
        fract_angle = len(u)/nframe
        if fract_angle<0.1:  ## ONLY USED very strong signals
            fract_angle=0.1
        fraction_angles.append(fract_angle)
    return (header_noFRAME,fraction_angles)

###############
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--rep1', type=str, default='',
                    help='INPUT FILE REP 1')
parser.add_argument('--rep2', type=str, default='',
                    help='INPUT FILE REP 2')
parser.add_argument('--start_res', type=int, default='UNKNOWN',
                    help='START residue')
parser.add_argument('--end_res', type=int, default='UNKNOWN',
                    help='END residue')
parser.add_argument('--cutoff', type=int, default='20',
                    help='Cutoff value to define BENDING. Default: 20 degrees')
parser.add_argument('--bcolor', type=str, default='bright',
                    help='Background color')
args = parser.parse_args()

ifile1 =  args.rep1
ifile2 =  args.rep2
background_color = args.bcolor
startRes = args.start_res
endRes = args.end_res
cutoff = args.cutoff

######################

data_rep1 = read_data(ifile1)
data_rep2 = read_data(ifile2)

# data_rep1 = data_rep1.head() ## FOR TESTING
# print(data_rep1)
del data_rep1['FRAME']

##############################################
# percent_a=[]
# percent_b=[]
# x_range=[]
# for i in np.arange(len(data_rep2)):
#     count_a=0
#     count_b=0

#     u1 = data_rep2.iloc[i]
#     # print(u1)
#     for j in np.arange(len(u1)):
#         if j<=34:
#             # print(i,j)
#             # print(u1[j])
#             if u1[j]>20:
#                 count_a+=1
#         else:
#             if u1[j]>20:
#                 count_b+=1

#     x_range.append(i)
#     percent_a.append(count_a)
#     percent_b.append(count_b)
#     # print(i,count_a/34,count_b/34)

# plt.plot(x_range,percent_a,label="CHAIN A")
# plt.plot(x_range,percent_b,label="CHAIN B")
# plt.legend()
##############################################

data_rep1['SYSTEM']='REP1'
data_rep2['SYSTEM']='REP2'

headers_rep1,fraction_angle_rep1 = fraction_of_angles(data_rep1,cutoff)
headers_rep2,fraction_angle_rep2 = fraction_of_angles(data_rep2,cutoff)

# print(headers_rep1)
combined = pd.DataFrame({'x':headers_rep1,
                        'y1':fraction_angle_rep1,
                        'y2':fraction_angle_rep2})

g = sns.lineplot(data=combined,
            x='x',
            y='y1',
            label='REP1',
            marker='o',
            # markersize=10,
            )
g.axvline(x = headers_rep1[18], zorder=0,color = 'black', linestyle = '--',label="Region Line") 

g = sns.lineplot(data=combined,
            x='x',
            y='y2',
            label='REP2',
            marker='o',
            # markersize=10,
            )
g.axvline(x = headers_rep1[52], zorder=0,color = 'black', linestyle = '--',label="Region Line") 

# del data_rep1['SYSTEM']  # Remove 'FRAME' for heatmap() plot
# sns.heatmap(data_rep1,vmin=cutoff,vmax=90)


# ### MISCELLANEOUS ###
plt.xticks(rotation=45)
# plt.ylim([0,1.1])
plt.suptitle("%s"%(ifile1[:-3]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(20,3)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
# plt.savefig("KDE%s"%(ifile[:-3]))
plt.show()

