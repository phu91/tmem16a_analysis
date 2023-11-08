#!/bin/bash 

module load gromacs 
XTC=$1
SYSTEM=$2	

echo q | gmx make_ndx -f step5_input.pdb -o index_new.ndx

echo 3 3 | /mnt/home/ptang/ceph/5-TMEM16A/delta-EAVK/1-MD/tmem16a_analysis/g_correlation_mpi -f $XTC -o correl_$SYSTEM.dat -s step5_input.pdb -n index_new.ndx -linear
echo 3 3 | /mnt/home/ptang/ceph/5-TMEM16A/delta-EAVK/1-MD/tmem16a_analysis/g_correlation_mpi -f $XTC -o mutual_$SYSTEM.dat -s step5_input.pdb -n index_new.ndx -linear -mi 
