# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
from MDAnalysis.analysis import hole2
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


def hydration_number(selection,pore_size):
	pore_water_up_5A = u.select_atoms("resname TIP3 and name OH2 and cyzone %s 5 0 (%s)"%(pore_size,selection),updating=True)
	pore_water_down_5A = u.select_atoms("resname TIP3 and name OH2 and cyzone %s 0 -5 (%s)"%(pore_size,selection),updating=True)
	return (len(pore_water_up_5A),len(pore_water_down_5A))

def pore_size(pore_center,frame,chain):
	# try:
	# 	os.mkdir("./hole2_output_%s"%(chain))
	# except OSError as error:
	# 	print(error)
	os.makedirs("./hole2_output_%s"%(chain),exist_ok=True)
	working_dir ='hole2_output_%s'%(chain)
	sel = u.select_atoms("protein")

	sel.write(working_dir+"/tmp.pdb")
	profiles = hole2.hole(working_dir+"/tmp.pdb",
					executable='~/hole2/exe/hole',
					outfile=working_dir+'/hole%i.out'%(frame),
					sphpdb_file=working_dir+'/hole%i.sph'%(frame),
					vdwradii_file=None,
					random_seed=31415,
					keep_files=True,
					cpoint=pore_center,
	   				cvect=[0,0,1])
	# hole2.create_vmd_surface(filename='hole2_output/hole%i.vmd'%(ts.frame),
	#                      sphpdb='hole2_output/hole%i.sph'%(ts.frame),
	#                     #  sph_process='~/hole2/exe/sph_process',
	#                   	  sos_triangle='~/hole2/exe/sos_triangle')

	# rxn_coords = profiles[0].rxn_coord
	# pore_length = rxn_coords[-1] - rxn_coords[0]
	pore_radius = profiles[0].radius
	pore_list = []
	with open(working_dir+'/radius_PROFILE_%s.out'%(chain),"a") as file_out:
		with open(working_dir+'/hole%i.out'%(frame),"r") as file_one:
			patrn = "sampled"
			for line in file_one:
				if re.search(patrn,line):
					line = line.split()
					pore_list = np.append(pore_list,(line[0],line[1]))
					file_out.write("%s\t%s\t%s\n"%(frame,line[0],line[1]))
				file_out.flush()
		file_one.close()
	file_out.close()

				# if pore_center[2]-0.8<=float(line[:12])<=pore_center[2]+0.8:
					# pore_list = np.append(pore_list,float(line[1:2]))
				# print(line)
				# print(len(line))
	# os.remove("tmp.pdb")
	# os.remove("hole1.*")
	# print(pore_list)
	# return pore_list

def pca_calculation(observable,n_pca):
	pca = PCA(n_components=n_pca)
	observable = StandardScaler().fit_transform(observable)
	principalComponents = pca.fit_transform(observable)
	principalDf = pd.DataFrame(data = principalComponents, columns = ['pc1', 'pc2'])
	return principalDf

def dssp_calculation(top_file,traj_file,frame):
	current_frame = md.load(traj_file,top=top_file)
	# print(current_frame)
	dssp1 = md.compute_dssp(current_frame)
	dssp1 = dssp1.tolist()
	H_content = []
	for i in range(len(dssp1)):
		n_residues = len(dssp1[i])
		count_H = dssp1[i].count('H')
		percent_H = count_H/n_residues
		H_content = np.append(H_content,percent_H)
	return H_content

top_file = sys.argv[1]
traj_file = sys.argv[2]

# print(top_file,traj_file)
chainA = "protein and segid PROA"
chainB = "protein and segid PROB"

# client = Client(n_workers=n_jobs)
u = mda.Universe(top_file,traj_file)


# helix6_628_A = "protein and resid 628 and segid PROA"
# helix6_655_A = "protein and resid 655 and segid PROA"
# helix6_668_A = "protein and resid 668 and segid PROA"
# helix6_full_A= "protein and resid 628 to 668 and segid PROA"
# helix6_628_B = "protein and resid 628 and segid PROB"
# helix6_655_B = "protein and resid 655 and segid PROB"
# helix6_668_B = "protein and resid 668 and segid PROB"
# helix6_full_B= "protein and resid 628 to 668 and segid PROB"
# helix6_full_A_sel = u.select_atoms(helix6_full_A)
# helix6_full_B_sel = u.select_atoms(helix6_full_B)

all = u.select_atoms("all")
chloride_a = u.select_atoms("(around 10 (resid 928 929 and resname CAL)) and (resname CLA)",updating=True)
chloride_b = u.select_atoms("(around 10 (resid 930 931 and resname CAL)) and (resname CLA)",updating=True)

with mda.Writer('cla_pore_a.dcd', all.n_atoms) as w:
    for ts in u.trajectory[::]:
        if len(chloride_a)!=0:
            print("Frame %s is collected."%(ts.frame))
            w.write(all)

with mda.Writer('cla_pore_b.dcd', all.n_atoms) as w:
    for ts in u.trajectory[::]:
        if len(chloride_b)!=0:
            print("Frame %s is collected."%(ts.frame))
            w.write(all)

# ofile = open("withCa-NoPIP2-PORE.dat",'w')
# ofile.write("frame PORE_A PORE_B\n")

# with open("withCa-NoPIP2-DSSP.dat","w+") as dssp_out:
# 	dssp_out.write("frame percentH_A percentH_B\n")
# 	for ts in u.trajectory[::]:
# 	    # print(ts.frame)
# 		helix6_full_A_sel.write("tmp_A_%s.pdb"%(ts.frame))
# 		dssp_result_A = dssp_calculation("tmp_A_%s.pdb"%(ts.frame),"tmp_A_%s.pdb"%(ts.frame),ts.frame)
# 		helix6_full_B_sel.write("tmp_B_%s.pdb"%(ts.frame))
# 		dssp_result_B = dssp_calculation("tmp_B_%s.pdb"%(ts.frame),"tmp_B_%s.pdb"%(ts.frame),ts.frame)

# 		dssp_out.write("%.3f\t%.3f\t%.3f\n"%(ts.frame,*dssp_result_A,*dssp_result_B))
# 		dssp_out.flush()
# 		os.remove("tmp_A_%s.pdb"%(ts.frame))
# 		os.remove("tmp_B_%s.pdb"%(ts.frame))

# dssp_out.close()

# 	pore_cen_A = u.select_atoms("%s"%(chainA)).center_of_geometry()
# 	pore_cen_B = u.select_atoms("%s"%(chainB)).center_of_geometry()
# 	pore_radius_A = pore_size(pore_cen_A,ts.frame,"A")
# 	pore_radius_B = pore_size(pore_cen_B,ts.frame,"B")


#	pore_water_A = hydration_number(chainA,pore_radius_avg_A)
#	pore_water_B = hydration_number(chainB,pore_radius_avg_B)

	# ofile.write("%i\t%.3f\t%.3f\n"%(ts.frame,pore_radius_avg_A,pore_radius_avg_B))     #,pore_water_A[0],pore_water_A[1],pore_water_B[0],pore_water_B[1]))
	# ofile.flush()
# ofile.close()
