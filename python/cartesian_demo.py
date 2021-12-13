import ismrmrd as mrd
import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt
from gadgetron.util.cfft import cifftn

# load data and mrd-header
dataset = mrd.Dataset('../data/converted_siemens_data.h5')
#dataset = mrd.Dataset('../data/ge_data/converted_data.h5')
exp_header = mrd.xsd.CreateFromDocument(dataset.read_xml_header())

# grab encoding information structure, and extract experiment dimensions
encoding_info = exp_header.encoding[0]

e_1 = encoding_info.encodedSpace.matrixSize.x
e_2 = encoding_info.encodedSpace.matrixSize.y # e_2 = encoding_info.encodingLimits.kspace_encoding_step_1.maximum
e_3 = encoding_info.encodedSpace.matrixSize.z 
e_4 = encoding_info.encodingLimits.average.maximum + 1
e_5 = encoding_info.encodingLimits.slice.maximum + 1
e_6 = encoding_info.encodingLimits.contrast.maximum + 1
e_7 = encoding_info.encodingLimits.phase.maximum + 1
e_8 = encoding_info.encodingLimits.repetition.maximum + 1
e_9 = encoding_info.encodingLimits.set.maximum + 1

# receiver channel information is in the acquisitionSystemInformation
num_channels = exp_header.acquisitionSystemInformation.receiverChannels

# total number of readouts in the experiment
num_acqs = dataset.number_of_acquisitions()

# print params to screen
exp_printout = [["samples", e_1],
                ["pe1", e_2],
                ["pe2", e_3],
                ["averages", e_4],
                ["slices", e_5],
                ["contrasts", e_6],
                ["phases", e_7],
                ["repetitions", e_8],
                ["sets", e_9],
                ["channels", num_channels],
                ["acqs", num_acqs]]
print(tabulate(exp_printout, headers=["experiment params", "values"]))

# create k-space array (channels first)
kspace = np.zeros((num_channels,e_1, e_2, e_3, e_4, e_5, e_6, e_7, e_8, e_9),dtype=np.complex64)

# line-by-line extraction and assignment to k-space
for ii in range(num_acqs):
    acq = dataset.read_acquisition(ii)
    if (ii==0):
        print(np.shape(acq.data)) # channel dimension is first

    # determine acquisition line in experiment
    kspace[:,:,
            acq.idx.kspace_encode_step_1,
            acq.idx.kspace_encode_step_2,
            acq.idx.average,
            acq.idx.slice,
            acq.idx.contrast,
            acq.idx.phase,
            acq.idx.repetition,
            acq.idx.set] = acq.data

# At this point, the data is ready for any python reconstruction!

# example 2D Cartesian reconstruction 
coil_images = cifftn(kspace, axes=[1,2])
image = np.squeeze(np.sqrt(np.sum(np.square(np.abs(coil_images[:,:,:,0,0,0,0,0,0,0])), axis=0)))
plt.imshow(image)
plt.show()

    
    
    
    


