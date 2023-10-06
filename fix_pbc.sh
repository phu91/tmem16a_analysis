#!/bin/sh

# INSTRUCTION
# This script will need these files:
	# 1. PDB
	# 2. PSF
	# 3. INDEX
	# 1. TPR
	# 1. XTC

############################################
# rechain.tlc 
############################################

#mol new step5_input.pdb
#mol addfile step5_input.psf

#set all [atomselect top all]

#$all move [trans x -90]
#$all move [trans y 180]

#set chaina [atomselect top "segname PROA"]
#set chainb [atomselect top "segname PROB"]
#set chainc [atomselect top "segname HETA"]
#set chaind [atomselect top "segname HETB"]

#$chaina set chain A
#$chainb set chain B
#$chainc set chain C
#$chaind set chain D

#$all writepdb NEWCHAIN.pdb
#$all writepsf NEWCHAIN.psf

#quit
############################################

### RUNNING COMMAND ###`
#./fix_pbc.sh TPR_FILE PRODUCTION_FILE SYSTEM_NAME
#####

 
TPR=$1
PRODUCTION=$2
SYSTEM=$3

module load modules/2.2-20230808 vmd

vmd -dispdev text -e rechain.tcl 

module unload modules/2.2-20230808 vmd

module load modules
module load gromacs

#gmx make_ndx -f $PDB -n index.ndx -o PBC_index.ndx << EOF
#3 | r SOD | r CAL | r CLA
#name 5 DRY 
#q
#EOF

gmx editconf -f NEWCHAIN.pdb -center 0 0 0 -o CENTERED_${SYSTEM}.pdb
echo 4 | gmx trjconv -f $PRODUCTION -s $TPR -pbc mol -n index.ndx -o tmp1.xtc
echo 0 4 | gmx trjconv -f tmp1.xtc -s CENTERED_${SYSTEM}.pdb -fit rot+trans -n index.ndx -o tmp2.xtc

mv tmp2.xtc PBC_${SYSTEM}.xtc
mv step5_input.psf PBC_${SYSTEM}.psf
rm tmp1.xtc

