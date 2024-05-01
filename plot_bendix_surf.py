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
    res1_start=res1[0]
    res1_end =res1[1]
    res2_start=res2[0]
    res2_end =res2[1]
    # print(res1_start)
    # print(res1_end)

    nres1=res1_end-res1_start+1
    nres2=res2_end-res2_start+1

    headers1=[]
    headers2=[]
    headers3=[]
    headers4=[]

    for i in np.arange(0,nres1*2):  # FOR two chain A and B
        if i <nres1:
            i = str(i+res1_start)+'A'
            headers1.append(i)
        else:
            i = str(i+res1_start-nres1)+'B'
            headers2.append(i)
    for i in np.arange(0,nres2*2):  # FOR two chain A and B
        if i <nres2:
            i = str(i+res2_start)+'A'
            headers3.append(i)
        else:
            i = str(i+res2_start-nres2)+'B'
            headers4.append(i)
    # headers=['FRAME']
    headers = headers1+headers3+headers2+headers4
    headers.insert(0,'FRAME')
    headers.insert(1,'FILE')
    # print(headers)
    return headers

def read_data(inputfile):
    df = pd.read_csv(inputfile,delim_whitespace=True,skiprows=14,
                    header=None,names=find_headers())
    df = df.drop_duplicates()
    # print(df)
    df['FILE']=inputfile
    return df

# def read_data(inputfile):
#     df = pd.read_csv(inputfile,delim_whitespace=True,skiprows=14,
#                     header=None)
#     df = df.drop_duplicates()
#     # print(df)
#     return df

def fraction_of_angles(dataframe,cutoff):
    print(dataframe)
    nframe = len(dataframe)
    header_noFRAME = find_headers()
    header_noFRAME.remove('FRAME')
    header_noFRAME.remove('0A')
    header_noFRAME.remove('0B')
    header_noFRAME.remove('FILE')
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
parser.add_argument('--input', nargs='+', help='<Required> Set flag', required=True)

parser.add_argument('--res1', type=int, default=[0,0], nargs='+',
                    help='START END Resid 1')

parser.add_argument('--res2', type=int, default=[0,0], nargs='+',
                    help='START END Resid 2')

parser.add_argument('--cutoff', type=int, default='20',
                    help='Cutoff value to define BENDING. Default: 20 degrees')

parser.add_argument('--bcolor', type=str, default='bright',
                    help='Background color')
args = parser.parse_args()

ifile = args.input
background_color = args.bcolor
res1 = args.res1
res2 = args.res2
cutoff = args.cutoff

######################
data_count = int(len(ifile))
# print(data_count)
data_collector = pd.DataFrame()
for i in range(data_count):
    data = read_data(ifile[i])
    # print(data)
    # data = data.dropna(axis=1)
    data_collector = pd.concat([data_collector, data],ignore_index=True)
data_collector = data_collector.drop(['0A','0B'], axis=1).fillna(0)  ## REMOVE EXTRA DATA
# print(data_collector)

u_mean_dataset = data_collector.groupby(['FRAME']).mean()  ## MEAN OF DATA SETS
# u_standardDev_dataset = data_collector.groupby(['FRAME']).std()
# u_standardDev_dataset = u_standardDev_dataset.fillna(0)
u_mean_dataset = u_mean_dataset.reset_index()
# u_standardDev_dataset = u_standardDev_dataset.reset_index()
headers_plot,fraction_angle_plot = fraction_of_angles(u_mean_dataset,cutoff)
# print(u_standardDev_dataset)
# print(headers_plot,fraction_angle_plot)
data_dict = pd.DataFrame({'x':headers_plot,
                        'y':fraction_angle_plot,
                        # 'yerr':u_standardDev_dataset.values()  ## STANDARD DEVIATION FROM MULTIPLE DATASETS
                        })

data_plot = pd.DataFrame(data_dict)
print(data_plot)
x_labels = data_plot.x
g = sns.lineplot(data=data_plot[:33],
x='x',
y='y')

# g.fill_between(u_mean_dataset, y-error, y+error)

# # # ### MISCELLANEOUS ###
# plt.xticks(rotation=72)
# # plt.ylim([0,1.1])
# plt.suptitle("%s"%(ifile1[:-3]),va='top')
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(25,3)   ## Wide x Height
# # plt.locator_params(axis='both', nbins=5)
# plt.tight_layout()
# # plt.savefig("KDE%s"%(ifile[:-3]))
plt.show()

