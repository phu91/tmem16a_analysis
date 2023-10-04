# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
from MDAnalysis.analysis import hole2, dihedrals
import os 
import warnings
warnings.filterwarnings("ignore")
import re
import mdtraj as md
import sys
# import dask
# import dask.multiprocessing
# dask.config.set(scheduler='processes')
# from dask.distributed import Client
from MDAnalysis.coordinates.memory import MemoryReader
import matplotlib.pyplot as plt

def cal_dihedral_omega(frame,residueList):
    omegas = [res.omega_selection() for res in residueList.residues]
    dihs = dihedrals.Dihedral(omegas).run(start=frame,stop=frame+1)
    return dihs.results.angles

def cal_dihedral_phi(frame,residueList):
    phis = [res.phi_selection() for res in residueList.residues]
    dihs = dihedrals.Dihedral(phis).run(start=frame,stop=frame+1)
    return dihs.results.angles

def cal_dihedral_psi(frame,residueList):
    psis = [res.psi_selection() for res in residueList.residues]
    dihs = dihedrals.Dihedral(psis).run(start=frame,stop=frame+1)
    return dihs.results.angles

top_file = sys.argv[1]   #PSF
traj_file = sys.argv[2]

# print(top_file,traj_file)
chainA = "protein and segid PROA"
chainB = "protein and segid PROB"

u = mda.Universe(top_file,traj_file)

all = u.select_atoms("all")
chloride = u.select_atoms("index 148924")

## ONLY SELECT CHAIN A WITH THE BINDING Chloride 148924
asn646 = u.select_atoms("protein and segid PROA and resid 646") 
lys584 = u.select_atoms("protein and segid PROA and resid 584")
gln645 = u.select_atoms("protein and segid PROA and resid 645")
lys641 = u.select_atoms("protein and segid PROA and resid 641")
asp550 = u.select_atoms("protein and segid PROA and resid 550")

inner_vestibule_a = u.select_atoms("protein and segid PROA and resid 646 584 645 641 550")
inner_vestibule_b = u.select_atoms("protein and segid PROB and resid 646 584 645 641 550")

# print(inner_vestibule.residues)
with open("dihedral_profile.dat","w+") as of:
    of.write("#protein and segid PROA | PROB and resid 646 584 645 641 550 for OMEGA | PHI | PSI (ordered)\n")
    for ts in u.trajectory[::10]:
        print("FRAME %s\n"%(ts.frame))
        frame = ts.frame
        omega_a = cal_dihedral_omega(ts.frame,inner_vestibule_a)
        omega_b = cal_dihedral_omega(ts.frame,inner_vestibule_b)
        phi_a   = cal_dihedral_phi(ts.frame,inner_vestibule_a)
        phi_b   = cal_dihedral_phi(ts.frame,inner_vestibule_b)
        psi_a   = cal_dihedral_psi(ts.frame,inner_vestibule_a)
        psi_b   = cal_dihedral_psi(ts.frame,inner_vestibule_b)

        for i,j,k,l,m,n in zip(omega_a,omega_b,phi_a,phi_b,psi_a,psi_b):
            of.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(frame,
                                                                   i[0],i[1],i[2],i[3],i[4],
                                                                   j[0],j[1],j[2],j[3],j[4],
                                                                   k[0],k[1],k[2],k[3],k[4],
                                                                   l[0],l[1],l[2],l[3],l[4],
                                                                   m[0],m[1],m[2],m[3],m[4],
                                                                   n[0],n[1],n[2],n[3],n[4]))
            of.flush()


# print(dihs.results)

# labels = ["Res %s"%(n) for n in inner_vestibule.residues.resids]
# print(labels)

# for ang,label in zip(dihs.angles.T,labels):
    # plt.plot(ang, label=label)

# plt.legend()
# plt.show()
# for ts in u.trajectory[::]:




# for ts in u.trajectory[::]:
    # pore_cen_A = u.select_atoms("%s"%(chainA)).center_of_geometry()
    # pore_cen_B = u.select_atoms("%s"%(chainB)).center_of_geometry()
    # pore_radius_A = pore_size(pore_cen_A,ts.frame,"A")
    # pore_radius_B = pore_size(pore_cen_B,ts.frame,"B")


#	pore_water_A = hydration_number(chainA,pore_radius_avg_A)
#	pore_water_B = hydration_number(chainB,pore_radius_avg_B)

