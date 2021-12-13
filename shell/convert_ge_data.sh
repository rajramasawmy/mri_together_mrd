#!/bin/sh 

# convert a pfile
ge2ismrmrd -v ../data/ge_data/P20480_GRE.7 -o ../data/ge_data/converted_ge_data.h5

# convert a scanarchive
ge2ismrmrd -v ../data/ge_data/ScanArchive_GRE.h5 -o ../data/ge_data/converted_ge_scan_archive.h5
