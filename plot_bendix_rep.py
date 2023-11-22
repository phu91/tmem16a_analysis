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
<<<<<<< HEAD

    order=['0A0B','3A3B','3A1B','3A0B']
    rep_list=['1','2']
    ## MAKE BOX PLOT for AVERAGE
    if kind=='box':
        box_pairs = [((("0A0B","A"),("0A0B","B"))),
                    ((("3A3B","A"),("3A3B","B"))),
                    ((("3A1B","A"),("3A1B","B"))),
                    ((("3A0B","A"),("3A0B","B"))),
                    ]

        g = sns.violinplot(data=data,
                        x='SYSTEM',
                        y='ANGLE',
                        hue='CHAIN',
                        split=True,
                        order=order,
                        palette='Set2',
                        zorder=2
                        )
        g.set_ylabel("MAXIMUM BENDING ANGLE")

        test_results = add_stat_annotation(g,data=data, x='SYSTEM', y='ANGLE',order=order,hue='CHAIN',
                                        box_pairs=box_pairs,
                                        test='Mann-Whitney', text_format='star',
                                        loc='inside', verbose=2)
=======
    
    order=['0A0B','3A3B','3A1B','3A0B']
    ## MAKE BOX PLOT for AVERAGE
    if kind=='box':
        box_pairs = [((("0A0B","A"),("0A0B","B"))),
                     ((("3A3B","A"),("3A3B","B"))),
                     ((("3A1B","A"),("3A1B","B"))),
                     ((("3A0B","A"),("3A0B","B"))),
                     ]

        g = sns.violinplot(data=data,
                     x='SYSTEM',
                     y='ANGLE',
                     hue='CHAIN',
                     split=True,
                     order=order,
                    palette='Set2',
                    zorder=2
                     )
        # g.axhline(y = 20, zorder=0,color = 'black', linestyle = '--',label="Straight") 
        # g.axhline(y = 90, zorder=0,color = 'red', linestyle = '--',label="Broken") 
        g.set_ylabel("MAXIMUM BENDING ANGLE")

        # test_results = add_stat_annotation(g,data=data, x='SYSTEM', y='ANGLE',order=order,hue='CHAIN',
        #                                    box_pairs=box_pairs,
        #                                    test='Mann-Whitney', text_format='star',
        #                                    loc='inside', verbose=2)
>>>>>>> 42ea1f6bf96f5c3b0dcc790f7d524f7f69868f55
        plt.suptitle("BOX_%s"%(ifile[:-4]),va='top')
    else:
        ### MAKE LINEPLOT FOR TIME
        fig,axes = plt.subplots(2,4,sharex=False,sharey=True,layout="constrained")
        # axes = axes.flatten()
        # print(axes)
        rep_list=['1','2']
        order=['0A0B','3A3B','3A1B','3A0B']
        data['TIME (ns)']=data['FRAME']*2
        # print(data.query("SYSTEM=='0A0B' and REP==1"))
        for row in np.arange(len(rep_list)):
            for col in np.arange(len(order)):
                # print(order[col], rep_list[row])
            # print(ax)
                g = sns.lineplot(data=data.query("SYSTEM=='%s' and REP==%s"%(order[col],rep_list[row])),
                                x='TIME (ns)',
                                y='ANGLE',
                                hue='CHAIN',
                                # style='SYSTEM',
                                palette='Set2',
                                ax=axes[row][col])
                g.set_title(" %s | REP %s"%(order[col], rep_list[row]))
                g.set_ylabel("MAXIMUM BENDING ANGLE")
                # g.axhline(y = 20, zorder=0,color = 'white', linestyle = '--',label="Straight") 
                # g.axhline(y = 90, zorder=0,color = 'red', linestyle = '--',label="Broken")     
        plt.suptitle("LINE_%s"%(ifile[:-4]),va='top')

    # ### MISCELLANEOUS ###
    plt.rcParams['ps.useafm'] = True
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rcParams['pdf.fonttype'] = 42
    plt.gcf().set_size_inches(7.5,5)   ## Wide x Height
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
                   names=['FRAME','CHAIN','ANGLE','SYSTEM','REP'])
print(data.isnull().values.any())
plot_bendix(dataFrame=data,kind=plot_kind,background=background_color)


