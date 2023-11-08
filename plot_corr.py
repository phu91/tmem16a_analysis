import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc

sns.set_context("notebook")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to RMSF profile')

args = parser.parse_args()
ifile = args.input

data = []
f = open ( ifile , 'r')
lines = f.readline()    # SKIP THE FIRST LINE
meta_info = lines.split(' ')
print(len(rows))

rows = meta_info[0]
cols = meta_info[1]

print(len(rows))
# lines = f.readlines()
# for line in lines:
#     data.extend(line.split())

# data.pop(-1)    # REMOVE THE LAST ']'
# # print(len(data))
# data = np.array(data)
# data_float = data.astype(np.float)
# # print(data_float)
# matrix = np.reshape(data_float,(rows,cols))
# print(matrix)
# ### MISCELLANEOUS ###
# plt.suptitle("%s"%(ifile[:-4]),va='top')
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,4)   ## Wide x Height
# # plt.locator_params(axis='both', nbins=5)
# # plt.tight_layout()
# plt.savefig("KDE%s"%(ifile[:-3]))
# plt.show()
