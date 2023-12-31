# Updated on 11/29/2022
import MDAnalysis as mda
import numpy as np
import os 
import warnings
warnings.filterwarnings("ignore")
import re
import mdtraj as md
import sys,stat
import argparse
import random as rnd
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import hole2
from tqdm import tqdm


def pore_radius_calculations(selection_frame,chain,frame,center_of_pore):
    # working_dir_file ='hole2_output_%s'%(systemname)
    # os.chmod(working_dir,stat.S_IRWXU)
    profiles = hole2.hole(selection_frame,
                executable='~/hole2/exe/hole',
                outfile=working_dir+'/hole_chain_%s_frame_%i.out'%(chain,frame),
                sphpdb_file=working_dir+'/hole_chain_%s_frame_%i.sph'%(chain,frame),
                vdwradii_file=None,
                random_seed=rnd.randrange(100000),
                keep_files=True,
                cpoint=center_of_pore,
                end_radius=40,
                cvect=[0,0,1])
    os.chmod(working_dir+'/hole_chain_%s_frame_%i.sph'%(chain,frame),stat.S_IRWXU)

    profiles_list = []
    with open(working_dir+'/hole_chain_%s_frame_%i.out'%(chain,frame),"r") as ifile:
        patrn = "sampled|mid-"
        for line in ifile:
            if re.search(patrn,line):
                line = line.split()
                profiles_list.append((line[0],line[1]))
    # print(frame,chain,profiles_list)
    return profiles_list

def draw_surface(chain,frame):
    # hole2.create_vmd_surface(filename=working_dir+'/hole_chain_%s_frame_%i.vmd'%(chain,frame),
    #                  sphpdb=working_dir+'/hole_chain_%s_frame_%i.sph'%(chain,frame),
    #                  sph_process='~/hole2/exe/sph_process',
    #               	 sos_triangle='~/hole2/exe/sos_triangle'
    #               	  )
    vmdfile=working_dir+'/hole_chain_%s_frame_%i.vmd'%(chain,frame)
    sphpdb=working_dir+'/hole_chain_%s_frame_%i.sph'%(chain,frame)
    sosfile=working_dir+'/hole_chain_%s_frame_%i.sos'%(chain,frame)
    sph_process='~/hole2/exe/sph_process'
    sos_triangle='~/hole2/exe/sos_triangle'
    os.system("%s -sos %s %s"%(sph_process,sphpdb,sosfile))
    os.system("%s -s <%s> %s"%(sos_triangle,sosfile,vmdfile))
    os.system("rm %s"%(sosfile))

def pore_radius_per_frame(pore_center,frame,chain):
    # working_dir ='hole2_output_%s'%(systemname)
    # os.chmod(working_dir,stat.S_IRWXU)
    if chain=='A':
        sel_frame = working_dir+"/tmpA_frame_%s.pdb"%(frame)
        sel = u.select_atoms("segid PROA").write(sel_frame)
        pore_profile=pore_radius_calculations(selection_frame=sel_frame,chain=chain,frame=frame,center_of_pore=pore_center)
        # print(pore_profile)
    else:
        sel_frame = working_dir+"/tmpB_frame_%s.pdb"%(frame)
        sel = u.select_atoms("segid PROB").write(sel_frame)
        pore_profile=pore_radius_calculations(selection_frame=sel_frame,chain=chain,frame=frame,center_of_pore=pore_center)
    # print(frame,chain,pore_profile)
    return (pore_profile)



# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate  Default = 1 frames')

parser.add_argument('--cutoff', type=int, default='5',
                    help='Cut-off distance to find Cl-. Default is 5A')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='System Name')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
cutoff = args.cutoff
systemname = args.system

u = mda.Universe(top_file,traj_file)
# os.system("rm hole2_output_%s/*.out"%(systemname))

try:
    os.makedirs("./hole2_output_%s"%(systemname))
except OSError as error:
    print(error)

working_dir ='hole2_output_%s'%(systemname)
os.chmod(os.path.abspath(working_dir),stat.S_IRWXU)

with open('PORE_PROFILE_%s.dat'%(systemname),'w+') as file_out:
    file_out.write("#frame chain z_coord radius\n")
    for ts in tqdm(u.trajectory[::traj_skip]):
        # print(ts.frame)
        com_pore_a = u.select_atoms("resid 550 551 584 641 645 646 700 701 and segid PROA", updating=True).center_of_geometry()
        com_pore_b = u.select_atoms("resid 550 551 584 641 645 646 700 701 and segid PROB", updating=True).center_of_geometry()
        file_out.write("# COM A %s %s\n"%(ts.frame,com_pore_a[2]))
        file_out.write("# COM B %s %s\n"%(ts.frame,com_pore_b[2]))

        pore_radius_A = pore_radius_per_frame(com_pore_a,ts.frame,"A")
        # print(pore_radius_A)
        for a,b in pore_radius_A:
            file_out.write("%s\t%s\t%s\t%s\n"%(ts.frame,'A',a,b))
            file_out.flush()

        pore_radius_B = pore_radius_per_frame(com_pore_b,ts.frame,"B")
        for a,b in pore_radius_B:
            file_out.write("%s\t%s\t%s\t%s\n"%(ts.frame,'B',a,b))
            file_out.flush()
        del pore_radius_A
        del pore_radius_B
        draw_surface('A',ts.frame)
        draw_surface('B',ts.frame)
        sel_frame_A_str = working_dir+"/tmpA_frame_%s.pdb"%(ts.frame)
        sel_frame_B_str = working_dir+"/tmpB_frame_%s.pdb"%(ts.frame)
        sel_frame_A_vmd_str = working_dir+'/hole_chain_A_frame_%i.vmd'%(ts.frame)
        sel_frame_B_vmd_str = working_dir+'/hole_chain_B_frame_%i.vmd'%(ts.frame)

        os.system("grep ATOM %s > tmp_%s_%s.pdb"%(sel_frame_A_str,systemname,ts.frame))
        os.system("grep ATOM %s >> tmp_%s_%s.pdb"%(sel_frame_B_str,systemname,ts.frame))

        os.system("mv tmp_%s_%s.pdb %s/frame_%s.pdb"%(systemname,ts.frame,working_dir,ts.frame))

        os.system("sed 1,1d %s > tmp_%s_%s.vmd"%(sel_frame_A_vmd_str,systemname,ts.frame))
        os.system("sed 1,1d %s >> tmp_%s_%s.vmd"%(sel_frame_B_vmd_str,systemname,ts.frame))

        os.system("mv tmp_%s_%s.vmd %s/frame_%s.vmd"%(systemname,ts.frame,working_dir,ts.frame))

        os.system("rm %s %s %s %s"%(sel_frame_A_str,sel_frame_B_str,sel_frame_A_vmd_str,sel_frame_B_vmd_str))
