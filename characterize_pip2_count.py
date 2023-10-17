# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader
# from MDAnalysis.analysis.leaflet import LeafletFinder
import math, sys
from MDAnalysis.lib.distances import distance_array
import argparse
from tqdm import tqdm

### FUNCTION
def lipid_enrichment(helix,PIP2):
    count=0
    # helix_com = helix.center_of_mass()
    pip2_closed = u.select_atoms("sphzone %s ( protein and ( %s ) and (name PLPI and name P)"%(rcutoff,helix))
    print(pip2_closed)
    # distances_a=distance_array(reference=PIP2.positions, 
    #                        configuration=helix.positions, 
    #                        box=u.dimensions, 
    #                        result=None, 
    #                        backend='OpenMP'
    # )

    # data_shape = np.shape(distances_a)
    # for row in range(data_shape[0]):
    #     for col in range(data_shape[1]):
    #         if distances_a[row,col]<=rcutoff:
    #             count+=1
    #             break  ## break the loop when one distance is less than rcut. MUST HAVE!!
    # # print(count)
    # return count


# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--begin', type=int, default=1,
                    help='Starting Frame. Default = 1 FRAME 1')

parser.add_argument('--end', type=int, default=-1,
                    help='Ending Frame. Default = -1 ALL FRAME')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate for Radial Distribution Function calculations. Default = 1 frames')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Add a system name to output file')

parser.add_argument('--rcut', type=int, default=3.5,
                    help='R cutoff value. Default = 3.5 Angstrom')

args = parser.parse_args()

top_file =  args.top
traj_file = args.traj
traj_begin = args.begin
traj_end = args.end
traj_skip = args.skip
systemname = args.system
rcutoff = args.rcut

u = mda.Universe(top_file,traj_file)

chainA_com = u.select_atoms("segid PROA",updating=True).center_of_geometry()
chainB_com = u.select_atoms("segid PROB",updating=True).center_of_geometry()

# Define HELIX UPDATED BASED ON BENDIX
helix_01a = "segid PROA and resid 327 to 360"
helix_02a = "segid PROA and resid 399 to 439"
helix_03a = "segid PROA and resid 478 to 518"
helix_04a = "segid PROA and resid 534 to 562"
helix_05a = "segid PROA and resid 572 to 598"
helix_06a = "segid PROA and resid 629 to 664"
helix_07a = "segid PROA and resid 692 to 713"
helix_08a = "segid PROA and resid 718 to 740"
helix_09a = "segid PROA and resid 753 to 780"
helix_10a = "segid PROA and resid 854 to 927"

helix_01b = "segid PROB and resid 327 to 360"
helix_02b = "segid PROB and resid 399 to 439"
helix_03b = "segid PROB and resid 478 to 518"
helix_04b = "segid PROB and resid 534 to 562"
helix_05b = "segid PROB and resid 572 to 598"
helix_06b = "segid PROB and resid 629 to 664"
helix_07b = "segid PROB and resid 692 to 713"
helix_08b = "segid PROB and resid 718 to 740"
helix_09b = "segid PROB and resid 753 to 780"
helix_10b = "segid PROB and resid 854 to 927"

helix_01a = "segid PROA and resid 327 to 360"
helix_02a = "segid PROA and resid 399 to 439"
helix_03a = "segid PROA and resid 478 to 518"
helix_04a = "segid PROA and resid 534 to 562"
helix_05a = "segid PROA and resid 572 to 598"
helix_06a = "segid PROA and resid 629 to 664"
helix_07a = "segid PROA and resid 692 to 713"
helix_08a = "segid PROA and resid 718 to 740"
helix_09a = "segid PROA and resid 753 to 780"
helix_10a = "segid PROA and resid 854 to 927"

helix_01b = "segid PROB and resid 327 to 360"
helix_02b = "segid PROB and resid 399 to 439"
helix_03b = "segid PROB and resid 478 to 518"
helix_04b = "segid PROB and resid 534 to 562"
helix_05b = "segid PROB and resid 572 to 598"
helix_06b = "segid PROB and resid 629 to 664"
helix_07b = "segid PROB and resid 692 to 713"
helix_08b = "segid PROB and resid 718 to 740"
helix_09b = "segid PROB and resid 753 to 780"
helix_10b = "segid PROB and resid 854 to 927"

helixa_list = [helix_01a,
helix_02a,
helix_03a,
helix_04a,
helix_05a,
helix_06a,
helix_07a,
helix_08a,
helix_09a,
helix_10a]

helixb_list = [helix_01b,
helix_02b,
helix_03b,
helix_04b,
helix_05b,
helix_06b,
helix_07b,
helix_08b,
helix_09b,
helix_10b]

helixa_list_label = ['helix_01a',
'helix_02a',
'helix_03a',
'helix_04a',
'helix_05a',
'helix_06a',
'helix_07a',
'helix_08a',
'helix_09a',
'helix_10a']

helixb_list_label = ['helix_01b',
'helix_02b',
'helix_03b',
'helix_04b',
'helix_05b',
'helix_06b',
'helix_07b',
'helix_08b',
'helix_09b',
'helix_10b']

# pip2 = u.select_atoms("resname PLPI and name P",updating=True)

lipids = u.select_atoms("resname DPPC DSPC DUPC PLPC PLPI POPC SLPC SOPC and name P")
lipid_com = u.select_atoms("resname DPPC DSPC DUPC PLPC PLPI POPC SLPC SOPC and name P").center_of_geometry()

############ Identify LEAFLETS
### Assuming the MB is wrap and centered with upperleaflet positioned above the oigin
upperleaflet = []
lowerleaflet = []
for lipid in lipids:
    if lipid.position[1]>= lipid_com[1]:
        # print(lipid.position[2],lipid_com[2])
        upperleaflet.append(lipid)
    else:
        lowerleaflet.append(lipid)
print("\nThere are %s lipids in the upperleaflet and %s lipids in the lower leaflet.\nTotal %s lipids"%(len(upperleaflet),len(lowerleaflet),len(upperleaflet)+len(lowerleaflet)))

# with open("PIP2_COUNT_%s.dat"%(systemname),'w+') as of:
#     of.write("#frame chain helix npip2\n")
#     for ts in tqdm(u.trajectory[traj_begin:traj_end:traj_skip]):
#         for ind, (helixa,helixb) in enumerate(zip(helixa_list,helixb_list)):
#             counta = lipid_enrichment(helixa,pip2)
#             countb = lipid_enrichment(helixb,pip2)
#             of.write("%s\t%s\t%s\t%s\n"%(ts.frame,"A",ind+1,counta))
#             of.flush()
#             of.write("%s\t%s\t%s\t%s\n"%(ts.frame,"B",ind+1,countb))
#             of.flush()

