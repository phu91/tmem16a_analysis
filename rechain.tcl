mol new step5_input.pdb
mol addfile step5_input.psf 

set all [atomselect top all]

$all move [trans x -90]
$all move [trans y 180]

set chaina [atomselect top "segname PROA"]
set chainb [atomselect top "segname PROB"]
set chainc [atomselect top "segname HETA"]
set chaind [atomselect top "segname HETB"]

$chaina set chain A 
$chainb set chain B 
$chainc set chain C 
$chaind set chain D 

$all writepdb NEWCHAIN.pdb
$all writepsf NEWCHAIN.psf 

quit
