## USAGE
# vmd -dispdev text -e analysis_scripts/change_chain.tcl [PDBFILE]
set all [atomselect top all]
set chaina [atomselect top "segname PROA and not resid 405 425"]
set chainb [atomselect top "segname PROB"]

$chaina set chain A
$chainb set chain B

$all writepdb corrected_chain.pdb
quit
