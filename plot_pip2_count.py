import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from statannot import add_stat_annotation

sns.set_context("notebook")
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to RMSF profile')

args = parser.parse_args()

ifile =  args.input

data = pd.read_csv(ifile,
                 names=['Frame','Helix','Chain','Enrichment'],
                 comment='#',
                 delim_whitespace=True)


chain_list = ['A','B']
color_list = sns.color_palette('Set2')

##### TEST
# df = sns.load_dataset("tips")
# x = "day"
# y = "total_bill"
# order = ['Sun', 'Thur', 'Fri', 'Sat']
# ax = sns.boxplot(data=df, x=x, y=y, order=order)
# test_results = add_stat_annotation(ax, data=df, x=x, y=y, order=order,
#                                    box_pairs=[("Thur", "Fri"), ("Thur", "Sat"), ("Fri", "Sun")],
#                                    test='Mann-Whitney', text_format='star',
#                                    loc='inside', verbose=2)
##### TEST

fig,axes = plt.subplots(2,1,sharex=False,sharey=False,layout="constrained")
axes = axes.flatten()
g1 = sns.boxplot(data=data, 
                x="Helix", 
                y="Enrichment",
                palette='Set2',
                hue="Chain",
                ax=axes[0])
g1.set_ylabel("Enrichment compared to Bulk")
g1.axhline(y = 1, color = 'black', linestyle = '--',label="Bulk") 

x2 = "Chain"
y2 = "Enrichment"
order = ['A','B']
g2 = sns.boxplot(data=data,
            x='Chain',
            y='Enrichment',
            palette='Set2',
            # hue='Chain'
            ax=axes[1]
            )
g2.set_ylabel("Enrichment compared to Bulk")
test_results = add_stat_annotation(axes[1],data=data, x=x2, y=y2,order=order,
                                   box_pairs=[("A", "B")],
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)

# ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,10)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
plt.savefig("%s"%(ifile[:-3]))
plt.show()
