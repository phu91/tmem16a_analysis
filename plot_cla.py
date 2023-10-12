import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc

sns.set(style="ticks", context="talk")
# plt.style.use("dark_background")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to Rg profile')

parser.add_argument('--pore_def', type=str, default='',
                    help='INPUT to PORE TOPOLOGY')

args = parser.parse_args()

ifile =  args.input
pore_def = args.pore_def

data = pd.read_csv(ifile,comment='#',
                   delim_whitespace=True,
                   names=['frame','Index','Found@Chain','Z','Distance'])

data['ns'] = data['frame']*0.1

pore_top = pd.read_csv(pore_def,
                  comment='#',
                  delim_whitespace=True,
                  names=['f','ca_min','ca_max','inner_min','inner_max','neck_min','neck_max','outer_min','outer_max'])

pore_top_avg = pore_top.mean().drop(['f'])

point_color = sns.cubehelix_palette(10,start=.5, rot=20, as_cmap=True)


sns.scatterplot(data=data.query("Distance <30"), 
            x="ns", 
            y="Z", 
            hue="Index",
            style="Found@Chain",
            alpha=0.7,
            legend='auto',
            palette=point_color,
            )

pore_cutoff = list(pore_top_avg.values)
# print(a)
pore_label = ["Ca2+ Binding","Inner Vestibule","Neck","Outer Vestibule"]
color_pore = sns.color_palette('Set2')

offset = np.linspace(0,100,5)
# print(offset)
for i,c in zip(np.arange(len(pore_label)),color_pore):
    plt.plot((offset[i],offset[i]),
             (pore_cutoff[2*i],pore_cutoff[2*i+1]),
             color=c,
             lw=10,
             label=pore_label[i],
             )


# ### MISCELLANEOUS ###
plt.ylim([-30,30])
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(10,7.5)
plt.locator_params(axis='y', nbins=6)
plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()
plt.savefig("%s"%(ifile[:-3]))
plt.show()
