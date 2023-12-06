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
    nres1=endRes1-startRes1+1
    nres2=endRes2-startRes2+1

    headers1=[]
    headers2=[]
    headers3=[]
    headers4=[]

    for i in np.arange(0,nres1*2):  # FOR two chain A and B
        if i <nres1:
            i = str(i+startRes1)+'A'
            headers1.append(i)
        else:
            i = str(i+startRes1-nres1)+'B'
            headers2.append(i)
    for i in np.arange(0,nres2*2):  # FOR two chain A and B
        if i <nres2:
            i = str(i+startRes2)+'A'
            headers3.append(i)
        else:
            i = str(i+startRes2-nres2)+'B'
            headers4.append(i)
    # headers=['FRAME']
    headers = headers1+headers3+headers2+headers4
    headers.insert(0,'FRAME')
    # print(headers)
    return headers

def read_data(inputfile):
    df = pd.read_csv(inputfile,delim_whitespace=True,skiprows=14,
                    header=None,names=find_headers())
    df = df.drop_duplicates()
    # print(df)
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
parser.add_argument('--rep3', type=str, default='',
                    help='INPUT FILE REP 3')
parser.add_argument('--rep4', type=str, default='',
                    help='INPUT FILE REP 4')
parser.add_argument('--rep5', type=str, default='',
                    help='INPUT FILE REP 5')
parser.add_argument('--start_res1', type=int, default='UNKNOWN',
                    help='START residue')
parser.add_argument('--end_res1', type=int, default='UNKNOWN',
                    help='END residue')
parser.add_argument('--start_res2', type=int, default='UNKNOWN',
                    help='START residue')
parser.add_argument('--end_res2', type=int, default='UNKNOWN',
                    help='END residue')
parser.add_argument('--cutoff', type=int, default='20',
                    help='Cutoff value to define BENDING. Default: 20 degrees')
parser.add_argument('--bcolor', type=str, default='bright',
                    help='Background color')
args = parser.parse_args()

ifile1 =  args.rep1
ifile2 =  args.rep2
ifile3 =  args.rep3
ifile4 =  args.rep4
ifile5 =  args.rep5
background_color = args.bcolor
startRes1 = args.start_res1
endRes1 = args.end_res1
startRes2 = args.start_res2
endRes2 = args.end_res2
cutoff = args.cutoff

######################

data_rep1 = read_data(ifile1)
data_rep2 = read_data(ifile2)
data_rep3 = read_data(ifile3)
data_rep4 = read_data(ifile4)
data_rep5 = read_data(ifile5)

# data_rep1 = data_rep1.head() ## FOR TESTING
# print(data_rep1)
# print(data_rep2)

del data_rep1['FRAME']
del data_rep2['FRAME']
del data_rep3['FRAME']
del data_rep4['FRAME']
del data_rep5['FRAME']


data_rep1['SYSTEM']='REP1'
data_rep2['SYSTEM']='REP2'
data_rep3['SYSTEM']='REP3'
data_rep4['SYSTEM']='REP4'
data_rep5['SYSTEM']='REP5'

headers_rep1,fraction_angle_rep1 = fraction_of_angles(data_rep1,cutoff)
# headers_rep2,fraction_angle_rep2 = fraction_of_angles(data_rep2,cutoff)
headers_rep3,fraction_angle_rep3 = fraction_of_angles(data_rep3,cutoff)
headers_rep4,fraction_angle_rep4 = fraction_of_angles(data_rep4,cutoff)
headers_rep5,fraction_angle_rep5 = fraction_of_angles(data_rep5,cutoff)


# # print(headers_rep1)
combined = pd.DataFrame({'x':headers_rep1,
                        'y1':fraction_angle_rep1,
                        # 'y2':fraction_angle_rep2,
                        'y3':fraction_angle_rep3,
                        'y4':fraction_angle_rep4,
                        'y5':fraction_angle_rep5,
                        })

print(combined)


# g = sns.lineplot(data=combined,
#             x='x',
#             y='y1',
#             label='REP1',
#             marker='o',
#             # markersize=10,
#             )
# g.axvline(x = headers_rep1[18], zorder=0,color = 'black', linestyle = '--',label="Region Line") 

# g = sns.lineplot(data=combined,
#             x='x',
#             y='y2',
#             label='REP2',
#             marker='o',
#             # markersize=10,
#             )
# g.axvline(x = headers_rep1[52], zorder=0,color = 'black', linestyle = '--',label="Region Line") 

# g = sns.lineplot(data=combined,
#             x='x',
#             y='y3',
#             label='REP3',
#             marker='o',
#             # markersize=10,
#             )
# g = sns.lineplot(data=combined,
#             x='x',
#             y='y4',
#             label='REP4',
#             marker='o',
#             # markersize=10,
#             )
# g = sns.lineplot(data=combined,
#             x='x',
#             y='y5',
#             label='REP5',
#             marker='o',
#             # markersize=10,
#             )


# del data_rep1['SYSTEM']  # Remove 'FRAME' for heatmap() plot
# g1=sns.heatmap(data_rep1,vmin=cutoff,vmax=90)
# del data_rep2['SYSTEM']  # Remove 'FRAME' for heatmap() plot
# g2=sns.heatmap(data_rep2,vmin=cutoff,vmax=90)
# del data_rep3['SYSTEM']  # Remove 'FRAME' for heatmap() plot
# g3=sns.heatmap(data_rep3,vmin=cutoff,vmax=90)
# del data_rep4['SYSTEM']  # Remove 'FRAME' for heatmap() plot
# g4=sns.heatmap(data_rep4,vmin=cutoff,vmax=90)
del data_rep5['SYSTEM']  # Remove 'FRAME' for heatmap() plot
g5=sns.heatmap(data_rep5,vmin=cutoff,vmax=90)

# # ### MISCELLANEOUS ###
plt.xticks(rotation=72)
# plt.ylim([0,1.1])
plt.suptitle("%s"%(ifile1[:-3]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(25,3)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
# plt.savefig("KDE%s"%(ifile[:-3]))
plt.show()

