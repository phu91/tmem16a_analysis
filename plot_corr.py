import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
# import plotly.express as px

sns.set_context("notebook")
# plt.style.use("dark_background")
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
df = df.fillna(0)
df = df.replace(2000,0)
print(df)

sns.heatmap(df)
plt.show()
