# mri_together_mrd
MRD data conversion tutorial @ MRI TOGETHER

https://mritogether.github.io/

Raj Ramasawmy @ NHLBI, NIH. 
Ahsan Javed @ NHLBI, NIH.
Adrienne Campbell-Washburn @ NHLBI, NIH.

# Contents 
- **data**
    - cartesian mrd
    - spiral mrd
    - siemens_raw_cartesian
    - ge_data_directory
    - bruker_data
- **python**
    - simple cartesian reconstruction demo
- **matlab**
    - spiral demo
- **shell**
    - siemens_raw_data_conversion
    - example folder conversion

# ISMRMRD Resources

## Source Code
https://github.com/ismrmrd/ismrmrd
https://github.com/ismrmrd/ge_to_ismrmrd
https://github.com/ismrmrd/siemens_to_ismrmrd

## ISMRMD API
https://ismrmrd.github.io/apidocs/1.5.0/

## Publication
https://onlinelibrary.wiley.com/doi/10.1002/mrm.26089 

# Previous Video Tutorials
Introduction to MRD/ISMRMRD

https://www.youtube.com/watch?v=LhyqkNXdsoc

Data Conversion and editing

https://www.youtube.com/watch?v=ggauvF13-OM&t=592s

# Installing
## Quick Install Guide Videos

[mrd-viewer](https://www.youtube.com/watch?v=7ocDY2qYQUA)

[ismrmrd & siemens converter](https://www.youtube.com/watch?v=rkPaLznuT0Q)

## Installation (everything)
- note that this is installing things to the user home folder - you might want to change this!
- This install creates a local directory for installed programs, which does not require sudo (but sudo is still required for package installations)
- You will have to add the ~/local directory to your path, i.e. something like..
```
export ISMRMRD_HOME=~/local
export LD_LIBRARY_PATH=${ISMRMRD_HOME}/lib:${LD_LIBRARY_PATH}
export PATH=~/local/bin:~/local/usr/bin:${ISMRMRD_HOME}/bin:${PATH}
```

```bash
cd ~
mkdir local

sudo apt-get install g++ cmake git
sudo apt-get install libboost-all-dev xsdcxx libxerces-c-dev libhdf5-serial-dev h5utils hdf5-tools libtinyxml-dev libxml2-dev libxslt1-dev
sudo apt-get -y install doxygen git-core graphviz libboost-all-dev libfftw3-dev libhdf5-serial-dev

pip3 install ismrmrd gadgetron numpy matplotlib tabulate
pip3 install git+https://github.com/ismrmrd/ismrmrdviewer.git

git clone https://github.com/ismrmrd/ismrmrd.git
mkdir ismrmrd/build && cd ismrmrd/build
cmake -DCMAKE_INSTALL_PREFIX=~/local/ ..
make install

cd ../../
git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git
mkdir siemens_to_ismrmrd/build && cd ismrmrd/build
cmake -DCMAKE_INSTALL_PREFIX=~/local/ ..
make install
```


# Siemens MRI Cartesian conversion and Python recon
```bash
cd shell/
sh convert_siemens_data.sh
cd ../python/
python3 cartesian_demo.py
```
# GE MRI Cartesian data conversion
To install ge2ismrmrd converter we need to install GE Orchestra on a linux distro with the correct libraries. We will be using a opensuse_15.2 image with orchestra installed for this demo. 
A more detailed installation with debuggin tips is available in the repo that is primarily maintained by Dr. Vinai Roopchansingh. (https://github.com/ismrmrd/ge_to_ismrmrd).

Please update the paths as needed in the install_ge_convertor.sh. We generally install in ~/local/ folder and thats what the install script would do as well.

```bash
cd shell/
sh install_ge_convertor.sh
source ~/.bashrc 
sh convert_ge_data.sh
```
# Siemens Spiral conversion (c++) and Bruker Cartesian (julia) Matlab reconstruction
See [matlab](matlab/README.md) folder for some more examples. 
