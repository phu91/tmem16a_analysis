import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from statannot import add_stat_annotation

sns.set_context("notebook")

def plot_bendix(dataFrame,kind,background):
    if background=='dark':
        plt.style.use("dark_background")

    order=['0A0B','3A3B','3A1B','3A0B']
    ## MAKE BOX PLOT for AVERAGE
    if kind=='box':
        box_pairs = [((("0A0B","A"),("0A0B","B"))),
                     ((("3A3B","A"),("3A3B","B"))),
                     ((("3A1B","A"),("3A1B","B"))),
                     ((("3A0B","A"),("3A0B","B"))),
                     ]

        g = sns.boxplot(data=data,
                     x='SYSTEM',
                     y='ANGLE',
                     hue='CHAIN',
                     order=order,
                    palette='Set2',
                    zorder=2
                     )
        # g.axhline(y = 20, zorder=0,color = 'black', linestyle = '--',label="Straight") 
        # g.axhline(y = 90, zorder=0,color = 'red', linestyle = '--',label="Broken") 
        g.set_ylabel("MAXIMUM BENDING ANGLE")

        test_results = add_stat_annotation(g,data=data, x='SYSTEM', y='ANGLE',order=order,hue='CHAIN',
                                           box_pairs=box_pairs,
                                           test='Mann-Whitney', text_format='star',
                                           loc='inside', verbose=2)
        plt.suptitle("BOX_%s"%(ifile[:-4]),va='top')
    else:
        ### MAKE LINEPLOT FOR TIME
        fig,axes = plt.subplots(2,2,sharex=True,sharey=True,layout="constrained")
        axes = axes.flatten()
        order=['0A0B','3A1B','3A3B','3A0B']
        data['TIME (ns)']=data['FRAME']*2
        for ind,ax in enumerate(axes):
            # print(ax)
            g = sns.lineplot(data=data.query("SYSTEM=='%s'"%(order[ind])),
                            x='TIME (ns)',
                            y='ANGLE',
                            hue='CHAIN',
                            # style='SYSTEM',
                            palette='Set2',
                            ax=ax)
            g.set_title("%s"%(order[ind]))
            g.set_ylabel("MAXIMUM BENDING ANGLE")
            # g.axhline(y = 20, zorder=0,color = 'white', linestyle = '--',label="Straight") 
            # g.axhline(y = 90, zorder=0,color = 'red', linestyle = '--',label="Broken")     
        plt.suptitle("LINE_%s"%(ifile[:-4]),va='top')

    # ### MISCELLANEOUS ###
    plt.rcParams['ps.useafm'] = True
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rcParams['pdf.fonttype'] = 42
    plt.gcf().set_size_inches(7.5,6)   ## Wide x Height
    # plt.locator_params(axis='both', nbins=5)
    plt.legend()
    # plt.tight_layout()
    if kind=='box':
        plt.savefig("DIHEDRAL_BOX_%seps"%(ifile[:-3]))
        plt.savefig("DIHEDRAL_BOX_%spng"%(ifile[:-3]))
    else:
        plt.savefig("DIHEDRAL_LINE_%seps"%(ifile[:-3]))        
        plt.savefig("DIHEDRAL_LINE_%spng"%(ifile[:-3]))
    plt.show()


parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT FILE')
parser.add_argument('--plot', type=str, default='box',
                    help='INPUT FILE')
parser.add_argument('--bcolor', type=str, default='bright',
                    help='INPUT FILE')
args = parser.parse_args()

ifile =  args.input
plot_kind = args.plot
background_color = args.bcolor

data = pd.read_csv(ifile,comment='#',
                   delim_whitespace=True,
                   names=['FRAME','CHAIN','ANGLE','SYSTEM'])

plot_bendix(dataFrame=data,kind=plot_kind,background=background_color)


