function [img_s] = mrd_recon_cartesian_bruker(dfile)
% function [data_struct] = recon_spiral(data_file,  noise_file)
%
% Simple spiral recon using MRD data format
%
% Requirements:
% ------------------------------
% The ISMRM sunrise toolbox
% ------------------------------
% from Michael Hansen @ NIH, USA
% http://hansenms.github.io/sunrise/
%
% available at
% https://github.com/hansenms/ismrm_sunrise_matlab
%
% with the following paths added:
% ismrm_sunrise_path = uigetdir(matlabroot, 'Select ISMRM sunrise dir');
%
%     addpath([ismrm_sunrise_path filesep 'ismrm_sunrise_matlab-master'])

% ------------------------------
% The ISMRMRD matlab class
% ------------------------------
% from the MRD/ISMRMRD community
%
% available at
% https://github.com/ismrmrd/ismrmrd/tree/master/matlab
% And to add:
% ismrmrd_path = uigetdir(matlabroot, 'Select ismrmrd-matlab dir');
%
%     addpath([ismrmrd_path])
%
%
% R Ramasawmy NHLBI June 2020

%%


%% Set up

make_nhlbi_utils;

% ============================================
% Load data
% ============================================

dsetin = ismrmrd.Dataset(dfile, 'dataset');
MRD_h = ismrmrd.xml.deserialize(dsetin.readxml);

if isfield(MRD_h.measurementInformation, 'protocolName')
    disp(['Reconstructing: ' MRD_h.measurementInformation.protocolName]);
end
raw_data    = dsetin.readAcquisition; % h5read(dfile,'/dataset/data');

% nhlbi_utils.plot_experiment(raw_data);

% ============================================
% Grab imaging parameters
% ============================================

% Ny  = MRD_h.encoding.reconSpace.matrixSize.y;
% Nz  = MRD_h.encoding.reconSpace.matrixSize.z;
% 
% averages        = 1 + (MRD_h.encoding.encodingLimits.average.maximum);
% slices          = 1 + (MRD_h.encoding.encodingLimits.slice.maximum);
% contrasts       = 1 + (MRD_h.encoding.encodingLimits.contrast.maximum);
% phases          = 1 + (MRD_h.encoding.encodingLimits.phase.maximum);
% sets            = 1 + (MRD_h.encoding.encodingLimits.set.maximum);
% reps            = 1 + (MRD_h.encoding.encodingLimits.repetition.maximum);

Ny  = MRD_h.encoding.encodedSpace.matrixSize.y;
Nz  = MRD_h.encoding.encodedSpace.matrixSize.z;
% 
averages        = 1 + single(max(raw_data.head.idx.average));
slices          = 1 + single(max(raw_data.head.idx.slice));
contrasts       = 1 + single(max(raw_data.head.idx.contrast));
phases          = 1 + single(max(raw_data.head.idx.phase));
sets            = 1 + single(max(raw_data.head.idx.set));
reps            = 1 + single(max(raw_data.head.idx.repetition));

samples         =      double(raw_data.head.number_of_samples(1)); % hoping the first one is a readout!
dt              =      raw_data.head.sample_time_us(1)*1e-6;
channels        =      double(raw_data.head.active_channels(1));

disp(' ');disp('### Experiment Dimensions ###');disp(' ');
Experiment_parameters = {'Samples', 'PE1', 'PE2', 'Averages', 'Slices', 'Contrasts', 'Phases', 'Repetitions', 'Sets', 'Channels'}';
Value = [samples Ny Nz averages slices contrasts phases reps sets channels]';
disp(table( Experiment_parameters,Value )); clear Experiment_parameters Value; disp(' ');

%% Noise checks

% if isempty(nfile)
dmtx = diag(ones(1,channels));
% else
%     dmtx = nhlbi_utils.noise_adjust(nfile, MRD_h, dt);
% end

%% Recon

nav_frames = raw_data.head.flagIsSet(raw_data.head.FLAGS.ACQ_IS_NAVIGATION_DATA);
% ref_frames = raw_data.head.flagIsSet(raw_data.head.FLAGS.ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING);
acq_frames = find(nav_frames==0);
ro_vec = 1:samples;

kspace = complex(zeros([MRD_h.encoding.encodedSpace.matrixSize.x Ny Nz averages slices contrasts phases reps sets channels],'single'));

for ii = acq_frames
    
    d1 = raw_data.data{ii};
    d1 = ismrm_apply_noise_decorrelation_mtx(d1, dmtx);
    
    kspace(ro_vec,...
        raw_data.head.idx.kspace_encode_step_1(ii)+1, ...
        raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
        raw_data.head.idx.average(ii)+1, ....
        raw_data.head.idx.slice(ii)+1 , ...
        raw_data.head.idx.contrast(ii)+1, ...
        raw_data.head.idx.phase(ii)+1, ...
        raw_data.head.idx.repetition(ii)+1, ...
        raw_data.head.idx.set(ii)+1, ...
        :) = d1;
    
end


%% transform to image-space

kspace = mean(kspace,4);

if Nz == 1
    % 2D
    cCoil_imgs = ismrm_transform_kspace_to_image(kspace, [1 2]);
else
    % 3D
    cCoil_imgs = ismrm_transform_kspace_to_image(kspace, [1 2 3]);
end

%% remove OS
% OverSampling = MRD_h.encoding.encodedSpace.fieldOfView_mm.x./MRD_h.encoding.reconSpace.fieldOfView_mm.x;
% if OverSampling == 2
%     
%     xdim = MRD_h.encoding.encodedSpace.matrixSize.x;
%     temp = reshape(cCoil_imgs, [xdim prod(size(cCoil_imgs))/xdim]);
%     temp = temp( (xdim/4): (xdim/4 + xdim/2 -1), : );
%     new_dims = size(cCoil_imgs); new_dims(1) = MRD_h.encoding.reconSpace.matrixSize.x;
%     cCoil_imgs = reshape(temp, new_dims);
%     
% end

%% Roemer coil combination "slice-by-slice"
figure,

dims = size(cCoil_imgs);
CCM_img =  zeros(dims(1:end-1));

for par = 1:Nz
    for slc = 1:slices
        for coc = 1:contrasts
            for phc = 1:phases
                for repc = 1:reps
                    for setc = 1:sets
                        
                        temp = cCoil_imgs(:,:,par,1,slc,coc,phc,repc,setc,:);
                        csm = ismrm_estimate_csm_walsh( squeeze( temp ) );
                        ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
                        CCM_img(:,:,par,1,slc,coc,phc,repc,setc) = abs( sum( squeeze( temp ) .* ccm_roemer_optimal, 3) );
                        
                        % ============================================
                        % View updated figure
                        % ============================================
                        
                        imshow(quick_scale(CCM_img(:,:,par,1,slc,coc,phc,repc,setc)),[0 4]);
                        drawstring = [];
                        if Nz > 1
                            drawstring = [drawstring 'Slice ' num2str(par) ' '];
                        end
                        if slices > 1
                            drawstring = [drawstring 'Slice ' num2str(slc) ' '];
                        end
                        if contrasts > 1
                            drawstring = [drawstring 'Contrast ' num2str(coc) ' '];
                        end
                        if phases > 1
                            drawstring = [drawstring 'Phase ' num2str(phc) ' '];
                        end
                        if reps > 1
                            drawstring = [drawstring 'Repetition ' num2str(repc) ' '];
                        end
                        if sets > 1
                            drawstring = [drawstring 'Set ' num2str(setc) ' '];
                        end
                        title(drawstring); drawnow;
                        
                    end
                end
            end
        end
    end
end

% Sum-of-Squares coil combine
% SOS_img = squeeze(sqrt(sum(cCoil_imgs.*conj(cCoil_imgs),length(size(cCoil_imgs)) )));
% CCM_img = squeeze(SOS_img);


%% Return variables
img = squeeze(CCM_img);

img_s.img = img;
img_s.header = MRD_h;


end


function y = quick_scale(x)
y = x./mean(x(:));
end