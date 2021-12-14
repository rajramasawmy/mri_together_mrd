# matlab demo

# Requirements: 
## The ISMRM sunrise / MIRT toolbox
from Michael Hansen's [course](http://hansenms.github.io/sunrise/)

available at https://github.com/hansenms/ismrm_sunrise_matlab

(which includes the [MIRT](https://web.eecs.umich.edu/~fessler/code/))

With the following paths added for Jeff Fessler's MIRT toolboxes as well:
```matlab
ismrm_sunrise_path = uigetdir(matlabroot, 'Select ISMRM sunrise dir');
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'vdspiral'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'utilities'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'graph'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'systems'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'wls'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'penalty'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'general'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'fbp'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'nufft'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'mri'])
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'mex' filesep 'v7']) 
addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master' filesep 'irt' filesep 'mri'])
```

## The ISMRMRD matlab class
From the MRD/ISMRMRD community

available at
https://github.com/ismrmrd/ismrmrd/tree/master/matlab

And to add:
```matlab
ismrmrd_path = uigetdir(matlabroot, 'Select ismrmrd-matlab dir');
addpath([ismrmrd_path])
```

# Process spiral data
## convert data

```bash
siemens_to_ismrmrd -f data/siemens_raw_spiral.dat -z 2 -o data/mrd_spiral.h5 
siemens_to_ismrmrd -f data/siemens_raw_spiral.dat -z 1 -o data/mrd_spiral_noise.h5 
```

## simple recon with matlab
```matlab
dirPath = [pwd filesep];
uncorrected_spiral_images = mrd_spiral_recon([dirPath mrd_spiral.h5],[dirPath mrd_spiral_noise.h5])
```
And using attached GIRF-corrected trajectories
```matlab
dirPath = [pwd filesep];
corrected_spiral_images = mrd_spiral_recon([dirPath mrd_spiral_traj.h5],[dirPath mrd_spiral_noise.h5])
```

# Cartesian Bruker conversion and recon

Bruker conversion using [mri-reco](https://magneticresonanceimaging.github.io/MRIReco.jl/latest/filehandling/#Conversion) in julia

Install Julia [for linux](https://julialang.org/downloads/platform/#linux_and_freebsd)

```bash
wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.0-linux-x86_64.tar.gz
tar zxvf julia-1.7.0-linux-x86_64.tar.gz
export PATH="$PATH:/path/to/<Julia directory>/bin"
```

Install mri-reco
https://magneticresonanceimaging.github.io/MRIReco.jl/latest/#Installation

```

```
