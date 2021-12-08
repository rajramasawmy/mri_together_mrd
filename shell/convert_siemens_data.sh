#!/bin/sh 

# convert the experiment
siemens_to_ismrmrd -f ../data/siemens_raw_cartesian.dat -z 2 -o ../data/converted_siemens_data.h5

# convert the included noise dependencies (for pre-whitening)
siemens_to_ismrmrd -f ../data/siemens_raw_cartesian.dat -z 2 -o ../data/converted_siemens_data.h5