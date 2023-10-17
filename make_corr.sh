#!/bin/sh

module load gromacs
vmd -e change_chain.tcl

gmx make_ndx -f corrected_chain.pdb -n index_2.ndx -o index_corr.ndx << EOF
chain A & a CA
chain B & a CA
q
EOF

echo 6 6 | ~/Downloads/g_correlation_1_0_3/g_correlation_src/g_correlation -f *_fixedPBC*.xtc -s corrected_chain.pdb -n index_corr.ndx -skip 10 -o correlA.dat

echo 7 7 | ~/Downloads/g_correlation_1_0_3/g_correlation_src/g_correlation -f *_fixedPBC*.xtc -s corrected_chain.pdb -n index_corr.ndx -skip 10 -o correlB.dat

