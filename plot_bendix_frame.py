import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("surf",delim_whitespace=True,skiprows=13,header=None)
df = df.drop_duplicates()

nres=34

INDEX = df[0]
COL = []

for col in df.columns:
    if col > nres:
        col = col-nres
    COL.append(col)

del df[0]
print(df)

g = sns.heatmap(df,vmin=25,vmax=90)
# g.xticklabels(COL)
plt.show()
