import MDAnalysis as mda
import pandas as pd
import numpy as np
import math, sys
import argparse
from MDAnalysis.coordinates.memory import MemoryReader
from tqdm import tqdm
from MDAnalysis.analysis import dihedrals
import warnings
# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

def cal_dihedral(residueList):
    results = []
    omegas = [res.omega_selection() for res in residueList.residues]
    omega = dihedrals.Dihedral(omegas).run(step=traj_skip)

    phis = [res.phi_selection() for res in residueList.residues]
    phi = dihedrals.Dihedral(phis).run(step=traj_skip)
    psis = [res.psi_selection() for res in residueList.residues]
    psi = dihedrals.Dihedral(psis).run(step=traj_skip)
    
    # print(np.shape(omega.results.angles))
    for ind,a in enumerate(residueList.residues.segids):
        for frame,(b,c,d) in enumerate(zip(omega.results.angles,phi.results.angles,psi.results.angles)):
            results.append((frame,a[-1],residueList.residues.resnames[ind],residueList.residues.resids[ind],b[ind],c[ind],d[ind],systemname))
    results_df = pd.DataFrame(results,columns=['#FRAME','CHAIN','RESIDUE NAME','RESID','OMEGA','PHI','PSI','SYSTEM'])
    return results_df

def extract_frames(startFrame,endFrame):
    whole_backbone = u.select_atoms("all")
    with mda.Writer("SEGMENT_FRAME_%s_TO_%s.pdb"%(startFrame,endFrame)) as pdb:
        pdb.write(whole_backbone)

    with mda.Writer("SEGMENT_FRAME_%s_TO_%s.xtc"%(startFrame,endFrame), all.n_atoms) as W:
        for ts in u.trajectory[startFrame:endFrame:traj_skip]:
            W.write(whole_backbone)
            
# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, PDB, etc.) file')

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

parser.add_argument('--npc', type=int, default=3,
                    help='Number of Principal Components. Default = 3')

parser.add_argument('--sel', type=str, default='name CA',
                    help='Using VMD selection with " ". Default: all')

parser.add_argument('--ref', type=str, default='step5_input.pdb',
                    help='If not provided, use default output (step5_input.pdb) from CHARMM-GUI as the reference (Default)')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
traj_begin = args.begin
traj_end = args.end
systemname = args.system
selectionatom=args.sel
ncomponent = args.npc
reference = args.ref

u = mda.Universe(top_file,traj_file,in_memory=True)
glu425 = u.select_atoms("protein and resid 425",updating=True)
asp879 = u.select_atoms("protein and resid 879",updating=True) 
asp884 = u.select_atoms("protein and resid 884",updating=True)

triad_list = [glu425,asp879,asp884]

if traj_end != -1:
    extract_frames(traj_begin,traj_end)
    u = mda.Universe("SEGMENT_FRAME_%s_TO_%s.pdb"%(traj_begin,traj_end),
                     "SEGMENT_FRAME_%s_TO_%s.xtc"%(traj_begin,traj_end),
                     in_memory=True)
    print("\n########################################################")
    print("TOP : SEGMENT_FRAME_%s_TO_%s.pdb"%(traj_begin,traj_end))
    print("TRAJ: SEGMENT_FRAME_%s_TO_%s.xtc"%(traj_begin,traj_end))
    # print("NUMBER OF ATOMS (ORIGINAL): %s"%(n_atom_origin))
    # print("NUMBER OF ATOMS (CURRENT) : %s"%(len(u.atoms)))
    print("NUMBER OF FRAMES: %s"%(len(u.trajectory)))
    print("########################################################\n")
else:
    print("\n########################################################")
    print("TOP : %s"%(top_file))
    print("TRAJ: %s"%(traj_file))
    # print("NUMBER OF ATOMS (ORIGINAL): %s"%(n_atom_origin))
    # print("NUMBER OF ATOMS (CURRENT) : %s"%(len(u.atoms)))
    print("NUMBER OF FRAMES: %s"%(len(u.trajectory)))
    print("########################################################\n")
    pass

dihedral_glu425 = cal_dihedral(glu425)
dihedral_asp879 = cal_dihedral(asp879)
dihedral_asp884 = cal_dihedral(asp884)

dihedral_profiles = pd.concat([dihedral_glu425,dihedral_asp879,dihedral_asp884],ignore_index=True)
dihedral_profiles.to_csv("DIHEDRAL_%s.csv"%(systemname),
                        sep='\t',
                        float_format='%.3f',
                        index=False
                        )


