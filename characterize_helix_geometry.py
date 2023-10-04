# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import math, sys
from MDAnalysis.analysis import helix_analysis as hel
import seaborn as sns 
import argparse

sns.set_context("paper")

# FUNCTIONS

def helenal_analysis(helix_name,skipping):
    h = hel.HELANAL(u,select=helix_name,ref_axis=[0, 0, 1])
    h.run(step=skipping)
    return h

# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate for Radial Distribution Function calculations. Default = 1 frames')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Add a system name to output file')

# parser.add_argument('--chain', type=str,
#                     help='Chain of protein')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
systemname = args.system
# chain_name = args.chain


u = mda.Universe(top_file,traj_file)

# Define HELIX
helix_01a = u.select_atoms("segid PROA and resid 327 to 360")
helix_02a = u.select_atoms("segid PROA and resid 399 to 439")
helix_03a = u.select_atoms("segid PROA and resid 493 to 518")
helix_04a = u.select_atoms("segid PROA and resid 534 to 562")
helix_05a = u.select_atoms("segid PROA and resid 572 to 598")
helix_06a = u.select_atoms("segid PROA and resid 629 to 664")
helix_07a = u.select_atoms("segid PROA and resid 692 to 713")
helix_08a = u.select_atoms("segid PROA and resid 714 to 738")
helix_09a = u.select_atoms("segid PROA and resid 753 to 775")
helix_10a = u.select_atoms("segid PROA and resid 854 to 881")

helix_01b = u.select_atoms("segid PROB and resid 327 to 360")
helix_02b = u.select_atoms("segid PROB and resid 399 to 439")
helix_03b = u.select_atoms("segid PROB and resid 493 to 518")
helix_04b = u.select_atoms("segid PROB and resid 534 to 562")
helix_05b = u.select_atoms("segid PROB and resid 572 to 598")
helix_06b = u.select_atoms("segid PROB and resid 629 to 664")
helix_07b = u.select_atoms("segid PROB and resid 692 to 713")
helix_08b = u.select_atoms("segid PROB and resid 714 to 738")
helix_09b = u.select_atoms("segid PROB and resid 753 to 775")
helix_10b = u.select_atoms("segid PROB and resid 854 to 881")

pip2 = u.select_atoms("resname PLPI24 and name P")

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

with open("RDF_profile_%s.dat"%(systemname),"w+") as rdf_out:
    rdf_out.write("#  helix chain radius gofR\n")
    for ind, (helixa,helixb) in enumerate(zip(helixa_list,helixb_list)):
        rdfa = rdf_pip2_analysis(helixa,traj_skip)
        rdfb = rdf_pip2_analysis(helixb,traj_skip)
        for bin,gofr in zip(rdfa.results.bins,rdfa.results.rdf):
            print("%s\t%s\t%s\t%s"%(ind+1,"A",bin,gofr))
            rdf_out.write("%s\t%s\t%s\t%s\n"%(ind+1,"A",bin,gofr))
            rdf_out.flush()
        for bin,gofr in zip(rdfb.results.bins,rdfb.results.rdf):
            print("%s\t%s\t%s\t%s"%(ind+1,"B",bin,gofr))
            rdf_out.write("%s\t%s\t%s\t%s\n"%(ind+1,"B",bin,gofr))
            rdf_out.flush()
