#!/bin/sh
pdb=$1
psf=$2
system=$3

foo1='$chaina'
foo2='$chainb'
foo3='$ca'

############################################
# rechain.tcl
############################################
echo "############################################"
echo "Create extract_CA.tcl file"
echo "############################################"

cat > extract_CA.tcl << EOF

mol new $pdb
mol addfile $psf

set chaina [atomselect top "segname PROA"]
set chainb [atomselect top "segname PROB"]

$foo1 set chain A
$foo2 set chain B

set ca [atomselect top "backbone"]
$foo3 writepdb backbone.pdb
$foo3 writepsf backbone.psf

quit
EOF

############################################
 
module load modules/2.2-20230808 vmd
vmd -dispdev text -e extract_CA.tcl 

############################################

module unload modules/2.2-20230808 vmd

module load modules
module load gromacs

# gmx editconf -f backbone.pdb -center 0 0 0 -o CENTERED_CA_${SYSTEM}.pdb
# echo 4 | gmx trjconv -f $PRODUCTION -s $TPR -pbc mol -n $INDEX -o tmp1.xtc
# echo 0 4 | gmx trjconv -f tmp1.xtc -s CENTERED_CA_${SYSTEM}.pdb -fit rot+trans -n $INDEX -o tmp2.xtc
# echo 1 4 | gmx trjconv -f tmp1.xtc -s CENTERED_CA_${SYSTEM}.pdb -fit rot+trans -n $INDEX -o tmp3.xtc

# mv tmp2.xtc PBC_CA_CENTERED_PROT_${SYSTEM}.xtc
# mv tmp3.xtc PBC_CA_CENTERED_MEMB_${SYSTEM}.xtc

# cp $PSF PBC_CA_${SYSTEM}.psf
# rm tmp1.xtc
# rm extract_CA.tcl
# rm NEWCHAIN.pdb NEWCHAIN.psf CENTERED_${SYSTEM}.pdb