# Updated on 11/29/2022
import pandas as pd
import numpy as np
import os 
import warnings
warnings.filterwarnings("ignore")
import re
import sys
import matplotlib.pyplot as plt
import seaborn as sns 

def get_histogram(filename,dataframe,edeges_list):
    with open(filename,"w+") as of:
        of.write("# edge1 edge2 half_edge mean std\n")
        for i in np.arange(nbin):
            u = dataframe.query("%s<=z<%s"%(edeges_list[i],edeges_list[i+1]))['r']
            of.write("%s\t%s\t%s\t%s\t%s\n"%(edeges_list[i],edeges_list[i+1],(edeges_list[i]+(edeges_list[i+1]-edeges_list[i])/2),u.mean(),u.std()))

data_file_1 = sys.argv[1]
data_file_2 = sys.argv[2]
nbin_arg = sys.argv[3]

df1 = pd.read_csv(data_file_1, delim_whitespace=True, names=['f','z','r'])
df2 = pd.read_csv(data_file_2, delim_whitespace=True, names=['f','z','r'])

# print(df.head())

zmin1 = df1['z'].min()
zmax1 = df1['z'].max()
zmin2 = df2['z'].min()
zmax2 = df2['z'].max()

nbin=int(nbin_arg)

edges1 = np.histogram_bin_edges(df1.z,bins=nbin,range=(zmin1,zmax1))
edges2 = np.histogram_bin_edges(df2.z,bins=nbin,range=(zmin2,zmax2))


get_histogram("radii_profile_histogram_A.dat",
              df1,
              edges1)

get_histogram("radii_profile_histogram_B.dat",
              df2,
              edges2)
