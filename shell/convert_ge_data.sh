#!/bin/sh 

# convert the experiment
ge2ismrmrd -v ../data/ge_data/P20480_GRE.7 -o ../data/ge_data/converted_ge_data.h5

# convert the included noise dependencies (for pre-whitening)
ge2ismrmrd -v ../data/ge_data/ScanArchive_GRE.h5 -o ../data/ge_data/converted_ge_scanarchine.h5
