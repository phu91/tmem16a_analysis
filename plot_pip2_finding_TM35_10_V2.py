import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
# import colorcet as cc

sns.set(style="ticks", context="notebook")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT profile')

args = parser.parse_args()

ifile =  args.input

data = pd.read_csv(ifile,comment='#',
delim_whitespace=True,names=['frame','chainA_887','chainB_887']
)

data['time']=data['frame']*0.2

chainA_pip2 = data.loc[data['chainA_887']==1]
chainB_pip2 = data.loc[data['chainB_887']==1]

print(len(chainA_pip2)/50000,len(chainB_pip2)/50000)