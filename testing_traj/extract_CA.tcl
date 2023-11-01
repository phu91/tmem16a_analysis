
mol new TEST.psf
mol addfile TEST.psf

set chaina [atomselect top "segname PROA"]
set chainb [atomselect top "segname PROB"]

$chaina set chain A
$chainb set chain B

set ca [atomselect top "backbone"]
$ca writepdb backbone.pdb
$ca writepsf backbone.psf

quit
