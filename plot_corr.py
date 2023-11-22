import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
# import plotly.express as px

sns.set_context("notebook")
plt.style.use("dark_background")
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to RMSF profile')

args = parser.parse_args()
ifile = args.input

with open(ifile ,'r') as f:
    lines = f.readline()    # SKIP THE FIRST LINE
    # print(lines)
    meta_info = lines.split(' ')

    rows = int(meta_info[0])
    cols = int(meta_info[2])
    data = []
    lines = f.readlines()
    for line in lines:
        data.extend(line.split())

    data.pop(-1)    # REMOVE THE LAST ']'
    # print(len(data))
    data = np.array(data,dtype=np.float64)
    matrix = np.reshape(data,(rows,cols))

resid_a = np.arange(73,928)
resid_b = np.arange(73,928)
resid_list = np.concatenate((resid_a,resid_b))
df = pd.DataFrame(matrix,columns=resid_list,index=resid_list)

# df_short = df.loc[224]
# df_short = df.loc[[73,100],:]
df_short = df
fig,ax = plt.subplots(1,1)
g = sns.heatmap(data=df_short,
                # vmin=0,
                # vmax=1,
                robust=True,
                cmap='RdPu_r')
g.set_xlabel("RESID")
g.set_ylabel("RESID")

# fig = px.imshow(matrix,x=resid_list,y=resid_list,text_auto=True)
# fig.show()
# g = ax.imshow(matrix)

# ax.set_xticks(np.arange(0,1710))
# ax.set_yticks(np.arange(0,1710))

# ax.set_xticklabels(resid_list,rotation=45)
# ax.set_yticklabels(resid_list,)
# fig.colorbar(g)


### MISCELLANEOUS ###
plt.locator_params(axis='y', nbins=30)
plt.locator_params(axis='x', nbins=30)
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,6)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
# plt.tight_layout()
plt.savefig("KDE%s"%(ifile[:-3]))
plt.show()
