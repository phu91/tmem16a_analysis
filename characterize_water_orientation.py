# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader
import math, sys
import argparse, math
from tqdm import tqdm
from MDAnalysis.analysis.waterdynamics import AngularDistribution as AD
from MDAnalysis.analysis.waterdynamics import WaterOrientationalRelaxation as WOR
from MDAnalysis.analysis.waterdynamics import SurvivalProbability as SP

# import multiprocessing
# from multiprocessing import Pool
# from functools import partial


# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate  Default = 1 frames')

parser.add_argument('--cutoff', type=float, default='3.5',
                    help='Cut-off distance to find water. Default is 3.5A')
                    
parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='System Name')

parser.add_argument('--rep', type=int, default=1,
                    help='Replica name. Default: 1')

args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
cutoff = args.cutoff
systemname = args.system
replica = args.rep

# u = mda.Universe(top_file,traj_file)

# water_chainA_ca_binding = u.select_atoms("name OH2 and (around %s (protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729 and backbone and segid PROA))"%(cutoff),updating=True)
# water_chainA_inner_vest = u.select_atoms("name OH2 and (around %s (protein and resid 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710 and backbone and segid PROA))"%(cutoff),updating=True)
# water_chainA_neck       = u.select_atoms("name OH2 and (around %s (protein and resid 505 506 507 508 540 542 543 590 591 592 593 635 636 711 712 and backbone and segid PROA))"%(cutoff),updating=True)
# water_chainA_outer_vest = u.select_atoms("name OH2 and (around %s (protein and resid 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716 and backbone and segid PROA))"%(cutoff),updating=True)
# water_chainB_ca_binding = u.select_atoms("name OH2 and (around %s (protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729 and backbone and segid PROB))"%(cutoff),updating=True)
# water_chainB_inner_vest = u.select_atoms("name OH2 and (around %s (protein and resid 544 545 546 547 548 549 585 586 587 588 589 637 638 640 642 644 645 705 706 707 708 709 710 and backbone and segid PROB))"%(cutoff),updating=True)
# water_chainB_neck       = u.select_atoms("name OH2 and (around %s (protein and resid 505 506 507 508 540 542 543 590 591 592 593 635 636 711 712 and backbone and segid PROB))"%(cutoff),updating=True)
# water_chainB_outer_vest = u.select_atoms("name OH2 and (around %s (protein and resid 509 510 512 513 514 518 520 530 534 535 539 594 595 616 618 619 630 631 632 633 634 716 and backbone and segid PROB))"%(cutoff),updating=True)


select = "name OH2 and (around %s (protein and resid 643 648 649 651 652 653 654 700 702 703 704 7025 727 729 and backbone and segid PROA)"%(cutoff)
universe = mda.Universe(top_file,traj_file)
sp = SP(universe, select, verbose=True)
sp.run(start=0, stop=6, tau_max=6)
tau_timeseries = sp.tau_timeseries
sp_timeseries = sp.sp_timeseries

# print in console
for tau, sp in zip(tau_timeseries, sp_timeseries):
      print("{time} {sp}".format(time=tau, sp=sp))

# plot
plt.xlabel('Time')
plt.ylabel('SP')
plt.title('Survival Probability')
plt.plot(tau_timeseries, sp_timeseries)
plt.show()