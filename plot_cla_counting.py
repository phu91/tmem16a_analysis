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

parser.add_argument('--cutoff', type=int, default='50',
                    help='Cutoff distance from the pore region to search for Cl-. Default 50 (A)')
args = parser.parse_args()

ifile =  args.input
cutoff = args.cutoff

data = pd.read_csv(ifile,comment='#',
delim_whitespace=True,
names=['FRAME','INDEX','CHAIN','ZPOS','POSITION_PORE','DCOM_PORE'])
data['ns']=data.FRAME*0.2
data['COUNT']=1
data_count = data.groupby(['INDEX','CHAIN','POSITION_PORE']).sum().reset_index()
data_count2 = data_count[['INDEX','CHAIN','POSITION_PORE','COUNT']].query("POSITION_PORE=='Inner_Vestibule' or POSITION_PORE=='Ca_Binding'")
data_count2['RETENTION']=data_count2['COUNT']*0.2

fig, axes = plt.subplots(2,1,sharey=True,sharex=True)
chain_list=['A','B']
color_list=['hot','cool']
chain_count=len(data_count2.CHAIN.unique())
if chain_count==1:
    if data_count2.CHAIN.unique()=='A':
        dummy_df =  {'INDEX':['None','None'],
            'CHAIN':'B', 
            'POSITION_PORE': ['Ca_Binding','Inner_Vestibule'],
            'COUNT': [0,0],
            'RETENTION':[0,0]}
        data_count2 = pd.concat([data_count2, pd.DataFrame(dummy_df)], ignore_index = True)
    else:
        dummy_df =  {'INDEX':['None','None'],
            'CHAIN':'A', 
            'POSITION_PORE': ['Ca_Binding','Inner_Vestibule'],
            'COUNT': [0,0],
            'RETENTION':[0,0]}
        data_count2 = pd.concat([data_count2, pd.DataFrame(dummy_df)], ignore_index = True)

# Draw a nested barplot by species and sex
for ind,ax in enumerate(axes):
    g = sns.barplot(
        data=data_count2.query("CHAIN=='%s'"%(chain_list[ind])), 
        x="POSITION_PORE", 
        y="RETENTION", hue="INDEX",
        palette="%s"%(color_list[ind]),        # alpha=.6,
        ax=ax
    )
    g.set_title("PORE CHAIN %s"%(chain_list[ind]),weight='bold')
    g.set_xlabel("")
    g.set_ylabel("Retention Time (ns)",weight='bold')
    g.set_ylabel(g.get_ylabel(),weight='bold')
    g.set_xticklabels(["Ca²⁺ Binding","Inner Vestibule"], minor=False,weight='bold')
    g.set_ylim([0,800])

    legend=ax.legend(loc="upper right",
    # bbox_to_anchor=(.5, -0.8), 
    ncol=2, 
    title="Chloride Index", 
    frameon=False,
    )
    for line in legend.get_lines():
        line.set_linewidth(0)
        line.set_marker('s')
        line.set_markersize(10)
    for text in legend.get_texts():
        text.set_fontweight("bold")
        text.set_fontsize(15)
    legend.get_title().set_fontweight("bold")

    for bar in g.patches:

      # Using Matplotlib's annotate function and
      # passing the coordinates where the annotation shall be done
      # x-coordinate: bar.get_x() + bar.get_width() / 2
      # y-coordinate: bar.get_height()
      # free space to be left to make graph pleasing: (0, 8)
      # ha and va stand for the horizontal and vertical alignment
        g.annotate(format(bar.get_height(), '.1f'), 
                        (bar.get_x() + bar.get_width() / 2, 
                        bar.get_height()), ha='center', va='center',
                        size=15, xytext=(0, 20),weight='bold',rotation=45,
                        textcoords='offset points')

# # # ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top',weight='bold')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
SCALE=1.5
plt.gcf().set_size_inches(7.5*SCALE,7.5*SCALE)
plt.tight_layout()
plt.savefig("%s_COUNT.png"%(ifile[:-4]))
# plt.show()
