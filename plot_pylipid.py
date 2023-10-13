import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("DATA.csv",
                 comment='#',
                 delim_whitespace=True,
                 names=['Residue',
                        'ResidueID',
                        'Duration',
                        'Durationstd',
                        'Occupancy',
                        'Occupancystd',
                        'LipidCount',
                        'LipidCountstd',
                        'BindingSiteID',
                        'BindingSiteDuration',
                        'BindingSiteOccupancy',
                        'BindingSiteLipidCount',
                        'BindingSiteKoff',
                        'BindingSiteKoffBootstrapavg',
                        'BindingSiteResidenceTime',
                        'BindingSiteRSquared',
                        'BindingSiteRSquaredBootstrapavg',
                        'BindingSiteSurfaceArea',
                        'BindingSitePoseRMSD']
)
df['system'] = 'NEW'


sns.lineplot(data=df,
             x='ResidueID',
             y='Duration',
             marker='o')

plt.show()