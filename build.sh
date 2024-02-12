#!/bin/sh
source "/Applications/ccp4-8.0/bin/ccp4.setup-sh"
export PATH=/usr/bin:$PATH
make -f Makefile.arch -j
./fix_library.sh

time ./sails \
 -pdbin data/5fji_dg.pdb \
 -mtzin data/5fji.mtz \
 -colin-fo FP,SIGFP \
 -colin-fc FWT,PHWT \
 -colin-free FREE \
 -predin data/a.map