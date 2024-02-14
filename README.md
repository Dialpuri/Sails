# Sails - Software for the Automated Identification of Linked Sugars

## Overview

Sails is a software package for automatic addition of glycans to protein models obtained through X-ray crystallography. 
Sails uses the output of GlycoFind to be able to build glycans quickly and accurately. 

## Authors

Jordan Dialpuri - York Structural Biology Laboratory, University of York

## Installation and setup
### Prerequisites

You must have
- Clipper
- MMDB2

installed, which come included in the CCP4. To ensure they are in your path, you must source the appropriate script. To do this run:

    source /opt/xtal/ccp4-X.X/bin/ccp4.setup-sh 
where X.X is your CCP4 version.


## Usage

To run the executable run:

    ./sails <PARAMS>

The available parameters are

    -pdbin <filename>
	-mtzin <filename>        
    -predin <filename>
    -colin-fo <colpath>
    -colin-hl <colpath> or -colin-phifom <colpath>
    -colin-fc <colpath>
    -colin-free <colpath>
    -resolution <float>

### Description 

pdbin - a model to add sugars to

mtzin - the input reflections data from the experiment e.g. 5fji.mtz

predin - the predicted map from GlycoFind

colin-fo - column headings for F<sub>obs</sub> e.g. FP, SIGFP

coln-hl - column headings for Hendrickson-Lattman coefficients (ABCD) e.g. sfcalc.A, sfcalc.B, sfcalc.C, sfcalc.D

colin-fc - column headings for F<sub>calc</sub> e.g FWT, PHWT

colin-free - column heading for Free R Flag

resolution - the resolution cutoff to use

## Development

### Apple Silicon

To build for Apple Silicon run

    ./build.sh

### Other

To build for other platformrs run

    make -j

