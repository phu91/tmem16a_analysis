import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc

sns.set(style="ticks", context="notebook")
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

### DEFINE PORE TOPOLOGY
PORE = data[['POSITION_PORE','ZPOS']]
Inner_Vestibule=PORE.loc[PORE.POSITION_PORE=="Inner_Vestibule"]
Ca_Binding=PORE.loc[PORE.POSITION_PORE=="Ca_Binding"]
Neck=PORE.loc[PORE.POSITION_PORE=="Neck"]
Outer_Vestibule=PORE.loc[PORE.POSITION_PORE=="Outer_Vestibule"]
Out_of_Zone=PORE.loc[PORE.POSITION_PORE=="Out_of_Zone"]

# print(Outer_Vestibule.min()[1])
pore_topology_color = sns.color_palette('Accent')
pore_topology_label = ['Inner Vestibule','Ca2+ Binding','Neck','Outer Vestibule']
pore_topology_name = [Inner_Vestibule,Ca_Binding,Neck,Outer_Vestibule]


for i,j,k in zip(pore_topology_name,pore_topology_label,pore_topology_color):
    plt.fill_between(data.query("DCOM_PORE<=%s"%(cutoff)).ns,i.min()[1],i.max()[1],alpha=1,color=k,label=j)
sns.scatterplot(data=data.query("DCOM_PORE<=%s"%(cutoff)), 
            x="ns", 
            y="ZPOS", 
            hue="INDEX",
            style="CHAIN",
            alpha=0.7,
            markers=True,
            legend='auto',
            palette='bright',
            linewidth=0.5,
            )
plt.xlabel("Time (ns)")
# # ### MISCELLANEOUS ###
# plt.ylim([-30,30])
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,5.5)
# plt.locator_params(axis='y', nbins=6)
plt.legend(loc='upper left',bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()
plt.savefig("%s.png"%(ifile[:-4]))
plt.show()
