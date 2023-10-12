#!/bin/sh

# INSTRUCTION
# This script will need these files:
	# 1. PDB
	# 2. PSF
	# 3. INDEX
	# 1. TPR
	# 1. XTC

### RUNNING COMMAND ###
#./fix_pbc.sh TPR PSF PDB INDEX PRODUCTION SYSTEM_NAME

TPR=$1
PSF=$2
PDB=$3
INDEX=$4
PRODUCTION=$5
SYSTEM=$6

foo1='$all'
foo2='$chaina'
foo3='$chainb'
foo4='$chainc'
foo5='$chaind'

############################################
# rechain.tcl
############################################
echo "############################################"
echo "Create rechain.tcl file"
echo "############################################"

cat > rechain.tcl << EOF
mol new $PDB
mol addfile $PSF

set all [atomselect top all]

$foo1 move [trans x -90]
$foo1 move [trans y 180]

set chaina [atomselect top "segname PROA"]
set chainb [atomselect top "segname PROB"]
set chainc [atomselect top "segname HETA"]
set chaind [atomselect top "segname HETB"]

$foo2 set chain A
$foo3 set chain B
$foo4 set chain C
$foo5 set chain D

$foo1 writepdb NEWCHAIN.pdb
$foo1 writepsf NEWCHAIN.psf

quit
EOF

############################################
 
module load modules/2.2-20230808 vmd
vmd -dispdev text -e rechain.tcl 

############################################

module unload modules/2.2-20230808 vmd

module load modules
module load gromacs

gmx editconf -f NEWCHAIN.pdb -center 0 0 0 -o CENTERED_${SYSTEM}.pdb
echo 4 | gmx trjconv -f $PRODUCTION -s $TPR -pbc mol -n $INDEX -o tmp1.xtc
echo 0 4 | gmx trjconv -f tmp1.xtc -s CENTERED_${SYSTEM}.pdb -fit rot+trans -n $INDEX -o tmp2.xtc

mv tmp2.xtc PBC_${SYSTEM}.xtc
mv $PSF PBC_${SYSTEM}.psf
rm tmp1.xtc
rm rechain.tcl
