#!/bin/csh -f

# echo "ok"
source /evgen/evgen_inclusive_e/setup
ln -s /evgen/evgen_inclusive_e/$1/dat 
/evgen/evgen_inclusive_e/$1/evgen_inclusive_e input.dat
