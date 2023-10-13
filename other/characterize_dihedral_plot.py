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


dihedral_data = sys.argv[1]

df = pd.read_csv(dihedral_data, 
                  delim_whitespace=True, 
                  comment='#',
                  names=['f','A_omega_646','A_omega_584','A_omega_645','A_omega_641','A_omega_550',
                         'B_omega_646','B_omega_584','B_omega_645','B_omega_641','B_omega_550',
                         'A_phi_646','A_phi_584','A_phi_645','A_phi_641','A_phi_550',
                         'B_phi_646','B_phi_584','B_phi_645','B_phi_641','B_phi_550',
                         'A_psi_646','A_psi_584','A_psi_645','A_psi_641','A_psi_550',
                         'B_psi_646','B_psi_584','B_psi_645','B_psi_641','B_psi_550'
                  ]
)

year_list=list(df.columns)
# print(year_list[1:])
df_long = pd.melt(df, value_vars=year_list[1:],value_name='test', ignore_index=False)
df_long['chain'] ='UK'
for ind,row in df_long.iterrows():
    if row['variable'].startswith('A_'):
        row['variable']=='A'
        print(row)
# print(df_long)

#protein and segid PROA | PROB and resid 646 584 645 641 550 for OMEGA | PHI | PSI (ordered)

# print(df.head())

# zmin1 = df1['z'].min()
# zmax1 = df1['z'].max()
# zmin2 = df2['z'].min()
# zmax2 = df2['z'].max()

# nbin=int(nbin_arg)

# edges1 = np.histogram_bin_edges(df1.z,bins=nbin,range=(zmin1,zmax1))
# edges2 = np.histogram_bin_edges(df2.z,bins=nbin,range=(zmin2,zmax2))


# get_histogram("radii_profile_histogram_1.dat",
#               df1,
#               edges1)

# get_histogram("radii_profile_histogram_2.dat",
#               df2,
#               edges2)
