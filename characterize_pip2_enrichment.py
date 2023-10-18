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
def lipid_enrichment(helix,pip2_str,rcutoff):
    count=0
    """
    Find PIP2 (resname PLPI* and P) within the cut-off distance and then find all lipids within the cut-off sphere. 
    Calculate the local density of each lipid in the cut-off
    Calculate entrichment = local density (PIP2) / bulk density (PIP2)
    """
    # selects all atoms a certain cutoff away from another selection,
    # e.g. around 3.5 protein selects all atoms not belonging to protein
    # that are within 3.5 Angstroms from the protein
    pip2_closed = u.select_atoms("(around %s (%s)) and (%s)"%(rcutoff,helix,pip2_str))
    lipid_closed = u.select_atoms("(around %s (%s)) and (name P)"%(rcutoff,helix))
    local_lipids = []
    for lipid in lipid_closed:
        if identify_lowerleaflet_lipids(lipid) is True:
            local_lipids.append(lipid)
    if len(local_lipids) !=0:
        local_density = len(pip2_closed)/len(local_lipids)
        enrichment = local_density/bulk_density()
        # print("Helix %s Enrichment %s"%(helix,enrichment))
        return enrichment
    else:
        # print("Use a larger R cut off!")
        return -1

def identify_lowerleaflet_lipids(lipidToCheck):
    """
    Identify LEAFLETS and calculate the bulk density of PIP2
    """
    lipid_com = u.select_atoms("resname DPPC DSPC DUPC PLPC PLPI* POPC SLPC SOPC and name P").center_of_geometry()
    if lipidToCheck.position[2]< lipid_com[1]:
        return True

def bulk_density():
    lowerlipid=[]
    total_pip2=[]
    lipid_com = u.select_atoms("resname DPPC DSPC DUPC PLPC PLPI* POPC SLPC SOPC and name P").center_of_geometry()
    lipids = u.select_atoms("resname DPPC DSPC DUPC PLPC PLPI* POPC SLPC SOPC and name P")
    # print(lipids.resnames)
    for lipid in lipids:
        if lipid.position[2]< lipid_com[1]:
            lowerlipid.append(lipid)
        if lipid.resname=="PLPI24":
            total_pip2.append(lipid)
    bulk_density_pip2 = len(total_pip2)/len(lowerlipid)
    # print("\n##### Total %s lipids and %s PIP2 in the lower leaflet. Bulk Density of PIP2 = %s"%(len(lowerlipid),len(total_pip2),bulk_density_pip2))
    return bulk_density_pip2

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

parser.add_argument('--rcutoff', type=float, default=3.5,
                    help='R cutoff value. Default = 3.5 Angstrom')

args = parser.parse_args()

top_file =  args.top
traj_file = args.traj
traj_begin = args.begin
traj_end = args.end
traj_skip = args.skip
systemname = args.system
rcutoff = args.rcutoff

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

pip2 = "resname PLPI* and name P"

if rcutoff is None:
    print("'\n##### A default RCUTOFF of 3.5 Angstrom is used\n")
else:
    print("\n##### RCUTOFF is %s Angstrom\n"%(rcutoff))

with open("PIP2_COUNT_%s.dat"%(systemname),'w+') as of:
    of.write("#frame helix chain enrichment\n")
    for ts in tqdm(u.trajectory[traj_begin:traj_end:traj_skip]):
        for ind, (helixa,helixb) in enumerate(zip(helixa_list,helixb_list)):
            lipid_enrichment_a = lipid_enrichment(helixa,pip2,rcutoff)
            lipid_enrichment_b = lipid_enrichment(helixb,pip2,rcutoff)
            of.write("%s\t%s\t%s\t%s\n"%(ts.frame,ind+1,"A",lipid_enrichment_a))
            of.flush()
            of.write("%s\t%s\t%s\t%s\n"%(ts.frame,ind+1,"B",lipid_enrichment_b))
            of.flush()

