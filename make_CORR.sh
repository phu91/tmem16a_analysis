#!/bin/sh

# INSTRUCTION
# This script will need these files:
	# 1. PDB
	# 2. PSF

### RUNNING COMMAND ###
#./make_corr.sh PSF PDB SYSTEM_NAME

PSF=$1
PDB=$2
TRJ=$3
SYSTEM=$4

foo1='$all'
foo2='$chaina'
foo3='$chainb'
foo4='$chainc'
foo5='$chaind'

module purge
module load modules/2.3-alpha4 vmd 

############################################
# rechain.tcl
############################################
echo "############################################"
echo "Create rechain.tcl file"
echo "############################################"

cat > rechain.tcl << EOF
mol new $PDB
mol addfile $PSF

set all [atomselect top "protein and name CA"]

$foo1 move [trans z 180]

set chaina [atomselect top "segname PROA and name CA"]
set chainb [atomselect top "segname PROB and name CA"]

$foo2 set chain A
$foo3 set chain B

$foo1 writepdb NEWCHAIN.pdb

quit
EOF

############################################

vmd -dispdev text -e rechain.tcl 

############################################
module purge 
module load modules/2.1.1-20230405 gromacs/singlegpunode-2022.4

echo q | gmx make_ndx -f NEWCHAIN.pdb -o index_corr.ndx

#gmx make_ndx -f NEWCHAIN.pdb -o index_corr.ndx <<EOF
#chain A & a CA
#chain B & a CA
#q
#EOF
############################################
sed -i -e 's/HSD/HIS/g' NEWCHAIN.pdb
echo 3 3 | /mnt/home/ptang/ceph/5-TMEM16A/delta-EAVK/1-MD/tmem16a_analysis/g_correlation_mpi -f $TRJ -o NETWORK_MutualInfo_$SYSTEM.dat -s $PDB -n index_corr.ndx -linear -mi
