import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
# import colorcet as cc

def find_chain_position(chain_com,pivot_com):
    position_X=0
    position_Y=0
    columns_names=pivot_com.columns
    # print(columns_names)
    for i in range(len(columns_names)):
        if columns_names[i]==chain_com.X[0]:
            position_X = i

    row_names=pivot_com.index
    print(row_names)
    for i in range(len(row_names)):
        if row_names[i]==chain_com.Y[0]:
            position_Y = i
    return position_X,position_Y

sns.set(style="ticks", context="notebook")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT profile')

args = parser.parse_args()

ifile =  args.input

data = pd.read_csv(ifile,comment='#',
                delim_whitespace=True,
                names=['FRAME','X','Y','Z_THICKNESS']
)
# print(data)
data_COM = pd.read_csv("%s_COM.dat"%(ifile[:-4]),comment='#',
                delim_whitespace=True,
                names=['FRAME','XA','YA','ZA','XB','YB','ZB']
)
# print(data_COM)
data = data.loc[data['FRAME']<=5]
data_average = data.groupby(['FRAME','X','Y']).mean().reset_index()
data_average = data_average[['X','Y','Z_THICKNESS']]


data_com_average = data_COM.mean()

chainA = data_com_average[['XA','YA']]
chainB = data_com_average[['XB','YB']]
chainA['Z_THICKNESS']=-1
chainB['Z_THICKNESS']=-1

chainA_reshaped = np.reshape(chainA,(1,3))
chainB_reshaped = np.reshape(chainB,(1,3))

chainA_reshaped_DF=pd.DataFrame(tuple(chainA_reshaped),columns=['X','Y','Z_THICKNESS'])
chainB_reshaped_DF=pd.DataFrame(tuple(chainB_reshaped),columns=['X','Y','Z_THICKNESS'])
# print(chainB_reshaped_DF)

data_com_average_empty_A = data_average.copy()
data_com_average_empty_A['Z_THICKNESS']=0
data_com_average_empty_A = pd.concat([chainA_reshaped_DF,data_com_average_empty_A]).sort_values(['X'])
# print(data_com_average_empty_A)
data_com_average_empty_B = data_average.copy()
data_com_average_empty_B['Z_THICKNESS']=0
data_com_average_empty_B = pd.concat([chainB_reshaped_DF,data_com_average_empty_B]).sort_values(['X'])

pivot_data_COM_A = pd.pivot_table(data_com_average_empty_A,columns='X',index='Y',values='Z_THICKNESS')
pivot_data_COM_B = pd.pivot_table(data_com_average_empty_B,columns='X',index='Y',values='Z_THICKNESS')
pivot_data = pd.pivot_table(data_average,columns='X',index='Y',values='Z_THICKNESS')

pos_com_X_A = find_chain_position(chainA_reshaped_DF,pivot_data_COM_A)[0]
pos_com_Y_A = find_chain_position(chainA_reshaped_DF,pivot_data_COM_A)[1]
pos_com_X_B = find_chain_position(chainB_reshaped_DF,pivot_data_COM_B)[0]
pos_com_Y_B = find_chain_position(chainB_reshaped_DF,pivot_data_COM_B)[1]
print(pos_com_X_A,pos_com_Y_A,pos_com_X_B,pos_com_Y_B)

im = plt.scatter(pos_com_X_A,pos_com_Y_A,s=300,label='CHAIN A')
im = plt.scatter(pos_com_X_B,pos_com_Y_B,s=300,label='CHAIN B')

im = plt.imshow(pivot_data,
            cmap='viridis',
            norm=None,
            aspect=None, 
            interpolation='spline36', 
            alpha=None, 
            vmin=0, vmax=None, 
            origin='lower', 
            extent=None)

cbar = plt.colorbar(im, extend='max', shrink=0.9,cmap='viridis')
cbar.set_label('Membrane Thickness')


plt.xlabel("X Coordinate",weight='bold')
plt.ylabel("Y Coordinate",weight='bold')
# Hide X and Y axes tick marks
plt.xticks([])
plt.yticks([])

# # plt.suptitle("%s"%(ifile))
# # plt.rcParams['ps.useafm'] = True
# # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# # plt.rcParams['pdf.fonttype'] = 42
# # plt.gcf().set_size_inches(7.5,5.5)
plt.legend(loc='upper right') #,bbox_to_anchor=(1.01, 1),borderaxespad=0)
# # plt.tight_layout()
# # # plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
plt.show()
