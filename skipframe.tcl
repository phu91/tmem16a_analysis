#!/bin/bash
XTC=$1
SKIP=$2

gmx trjconv -f $XTC -skip $SKIP -o SKIP$SKIP.xtc
