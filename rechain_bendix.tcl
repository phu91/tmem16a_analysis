mol new PBC_3A0B_R1.psf
mol addfile PBC_CENTERED_PROT_3A0B_R1.xtc step 10 waitfor all
set a [atomselect top "segname PROA"]
set b [atomselect top "segname PROB"]
$a set chain A
$b set chain B
