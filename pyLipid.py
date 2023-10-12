import pylipid
from pylipid.api import LipidInteraction
import argparse

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--top', type=str, default='',
                    help='TOPOLOGY FILE')

parser.add_argument('--traj', type=str, default='',
                    help='TRAJECTORY FILE')

parser.add_argument('--lipid', type=str, default='',
                    help='LIPID NAME')
args = parser.parse_args()
topfile  = args.top
trajfile = args.traj
lipidname= args.lipid

trajfile_list = [trajfile]
topfile_list = [topfile]
lipid = lipidname
stride = 10
cutoffs = [0.475, 0.8] # use of dual-cutoff
nprot = 1  # num. of protein copies in the system. if the simulation system has N copies of receptors,
           # "nprot=N" will report interactions averaged from the N copies, but "nprot=1"
           # will ask pylipid to report interaction data for each copy.
timeunit = 'ns' # micro-second
save_dir = None  # if None, pylipid data will be saved at current working directory.

# initialize
li = LipidInteraction(trajfile_list, topfile_list=topfile_list, cutoffs=cutoffs, lipid=lipid,
                      nprot=nprot, save_dir=save_dir,stride=stride)
                      
li.collect_residue_contacts()
# print(li.dataset)
durations_res = li.compute_residue_duration()
# print(len(li.residue_list), len(li.trajfile_list))
# print(len(durations), len(durations[0]))
# print(durations[0])
occupancies_res = li.compute_residue_occupancy()
lipidcounts_res = li.compute_residue_lipidcount()
# print(len(occupancies), len(occupancies[0]))
# print(len(lipidcounts), len(lipidcounts[0]))
# print(li.dataset)

# koffs, res_times = li.compute_residue_koff(plot_data=True, fig_close=False)
# fig_close=True will close all the figures. As it can generate substantial amount of figures
# (one for each residue), leaving all the figures open can consume considerable amount of memory.

# li.compute_residue_koff(residue_id=[10,15,25,35])
# li.compute_residue_koff(residue_id=10)
# li.dataset.to_csv("res.csv",sep='\t')

node_list, modularity = li.compute_binding_nodes(threshold=4, print_data=True)

# threshold=4 decides that binding sites should contain at least 4 residues. This is particularly
# helpful when itneraction samplings are not sufficient, in which case false positively correlation
# is likely to happen among 2 or 3 residues.
#
# print_data=True will print the residues for each binding site. It can be quite verbose.

durations_bs = li.compute_site_duration()
occupancies_bs = li.compute_site_occupancy()
lipidcounts_bs = li.compute_site_lipidcount()

koff_BS, res_time_BS = li.compute_site_koff(plot_data=True)

surface_area = li.compute_surface_area(plot_data=True)
surface_area.to_csv("SURFACE_AREA.csv",sep='\t')

pose_pool, pose_rmsd = li.analyze_bound_poses(binding_site_id=None, pose_format="pdb",
                                             n_top_poses=3, n_clusters='auto') # None means analysis for all residues
pose_rmsd.to_csv("POSE.csv")
li.dataset.to_csv("DATA.csv",sep='\t')

for item in ["Residence Time", "Duration", "Occupancy", "Lipid Count"]:
    li.plot(item=item) # plot values as a function of time

li.plot_logo(item="Residence Time")
li.write_site_info(sort_residue="Residence Time") # write binding site information in a txt file
li.save_coordinate(item="Residence Time") # write protein coordinate in pdb format with the b factor column
                                          # showing 'Residence Time' values