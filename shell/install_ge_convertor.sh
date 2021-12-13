#! /bin/bash

mkdir ~/mrprogs/
mkdir ~/local/
mkdir ~/local/ismrmrd
mkdir ~/local/ge-tools

#Define PATHS
#Set SDKTOP to GEorchestra, ISMRMRD and GETOOLs home to user local folders this 
export SDKTOP=/usr/local/devGE/Orchestra
export ISMRMRD_HOME=~/local/ismrmrd
export GE_TOOLS_HOME=~/local/ge-tools
export GADGETRON_HOME=~/local
export HDF5_ROOT=$SDKTOP/3p
export PATH="$ISMRMRD_HOME/bin/:$GE_TOOLS_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$ISMRMRD_HOME/lib/:$GE_TOOLS_HOME/lib/:$LD_LIBRARY_PATH"

echo 'export PATH="$ISMRMRD_HOME/bin/:$GE_TOOLS_HOME/bin:$PATH"' >> ~/.bashrc 
echo 'export LD_LIBRARY_PATH="$ISMRMRD_HOME/lib/:$GE_TOOLS_HOME/lib/:$LD_LIBRARY_PATH"' >> ~/.bashrc 
echo 'export ISMRMRD_HOME=~/local/ismrmrd'      >> ~/.bashrc
echo 'export GE_TOOLS_HOME=~/local/ge-tools'	>> ~/.bashrc
echo 'export GADGETRON_HOME=~/local'	        >> ~/.bashrc

source ~/.bashrc


mkdir ~/mrprogs/
mkdir -p $ISMRMRD_HOME
mkdir -p $GE_TOOLS_HOME

# Get ISMRMRD
cd ~/mrprogs/
git clone https://github.com/ismrmrd/ismrmrd

cd ismrmrd/
mkdir build
cd build/
cmake -D build4GE=ON -D CMAKE_INSTALL_PREFIX=$ISMRMRD_HOME ..
make install -j8
cd ../

cd ~/mrprogs/
git clone https://github.com/ismrmrd/ge_to_ismrmrd.git

cd ge_to_ismrmrd/
mkdir build
cd build/
cmake -D CMAKE_INSTALL_PREFIX=$GE_TOOLS_HOME ..
make install -j8
cd ../
