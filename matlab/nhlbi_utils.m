% NHLBI 
% Offline recon toolbox
% (This is real nasty "class"-ing.)

% -- Search for local functions --
function fh = nhlbi_utils
fh = localfunctions;
end

% -- edit call --
function edit_toolbox
edit(mfilename)
end

function [weights] = ismrm_calculate_weights(traj_info, method)
% function [weights] = ismrm_calculate_weights(traj_info, method)
%
% Assuming traj_info.traj is scaled to kr_max (not always needed)
% and has dims: 
%       traj_info.traj = [samples shots 2]
% ------------------------------
%
% method options / requirements
% ------------------------------
% method = 'Pipe'; <default> 
%       traj_info.traj
%       traj_info.FOV [X]
%       traj_info.res [X Y]
%        
% method = 'Meyer';
%       traj_info.traj
%       traj_info.grad
%
% method = 'Hoge';
%       traj_info.traj
%       traj_info.grad
%
% method = 'Voronoi';
%       traj_info.traj
%

if nargin < 1
    help ismrm_calculate_weights
    return
elseif nargin == 1
    method = 'Pipe';
end


%% UNPACK
    
traj = traj_info.traj;
[samples, shots, dims] = size(traj);
traj = reshape(traj, [samples*shots 2]);

if dims > 2
    warning('Only Pipe can deal with 3D');
end


%% METHODS
switch method
    case 'Pipe'
        % =============================
        % Pipe
        % =============================
        % 
        % JG Pipe & P Menon. Sampling Density Compensation in MRI:
        % Rationale and an Iterative Numerical Solution. MRM; 1999: 41:179?186
        %
        % Fesslers's IRT > Example from CWR implementation
        % requires a few inputs
        
        kr_max = 1./(2*traj_info.FOV/traj_info.res(1));
        traj = (traj./max(traj(:))).*kr_max;
        
        mask = true( traj_info.res );
        sizeMask = size(mask);
        nufft_args = {sizeMask, [6 6], 2*sizeMask, sizeMask/2, 'table', 2^12, 'minmax:kb'};
        G = Gmri(traj, mask, 'fov', traj_info.FOV, 'basis', {'dirac'}, 'nufft', nufft_args);
        
        weights = abs(mri_density_comp(traj, 'pipe','G',G.arg.Gnufft));
        
    case 'Hoge'
        % =============================
        % Hoge
        % =============================
        %
        % RD Hoge, et al. Density Compensation Functions for Spiral MRI. MRM; 1997:38-117-120 :: Eq [11]
        % https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.1910380117 
        % I made this up so.. 
        % Ramasawmy
        
        traj = reshape(traj, [samples*shots 2]);
        
        grad = traj_info.grad;
        grad = reshape(grad,[samples*shots 2]);
        
        grad = complex(grad(:,1),grad(:,2));
        traj = complex(traj(:,1),traj(:,2));
        weights = abs(traj(:)) .* abs(grad(:)) .* abs(cos( angle(grad(:)) - angle(traj(:)) )); 
        
    case 'Voronoi'
        % =============================
        % Voronoi
        % =============================
        %
        % Issues with overlapping values
        % Ramasawmy
       
        % interface
        k = complex(traj(:,1),traj(:,2));
        
        % function area = voronoidens(k);
        %
        % input:  k = kx + i ky is the  k-space trajectory
        % output: area of cells for each point
        %           (if point doesn't have neighbors the area is NaN)
        %
        % Written by John Pauly, modified by Michael Lustig
        
        
        r = max(abs(k(:)));
        k = [k(:); r*1.005*exp(1j*2*pi*(1:256).'/256)]; % radius to avoid blowing up
        kx = real(k);
        ky = imag(k);
        
        % uncomment these to plot voronoi diagram
        % [vx, vy] = voronoi(kx,ky);
        %figure, plot(kx,ky,'r.',vx,vy,'b-'); axis equal
        
        kxy = [kx(:),ky(:)];
        % returns vertices and cells of voronoi diagram
        [V,C] = voronoin(kxy);
        
        % weights = [];
        weights = zeros(length(kxy),1);
        
        for j = 1:length(kxy)
            if ~isempty(C{j})
                x = V(C{j},1); y = V(C{j},2); lxy = length(x);
                A = abs(sum( 0.5*(x([2:lxy 1]) - x(:)).*(y([2:lxy 1]) + y(:))));
            else
                A = inf;
            end
            % weights = [weights A];
            weights(j) = A;
        end
        
        weights = weights(1:end-256);
%         weights = weights(:)/sum(weights(:));
        
    case 'Meyer'
        % =============================
        % Meyer
        % =============================
        % spiral-out happy
        % Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
        
        grad = traj_info.grad;
        grad = reshape(grad,[samples*shots 2]);
        
        grad = complex(grad(:,1),grad(:,2));
        traj = complex(traj(:,1),traj(:,2));
        weights = abs(grad(:)) .* abs(sin( angle(grad(:)) - angle(traj(:)) )); 
    
    otherwise
        help ismrm_calculate_weights;
        return;
                
end

%% SCALE for NUFFT

% normalisation
weights = weights(:)/(sum(weights(:)));

% SNR-scaled?
% weights = weights*sqrt(prod(traj_info.res));

% reshape back?
% weights = reshape(weights, [samples shots]);

end

function plot_experiment(raw_data, subset_indices)

hl = [];
temp = fieldnames(raw_data.head.idx);
figure('Name', 'Experiment header');
for i = 1:9 % dont plot user vector
    h.(matlab.lang.makeValidName(['x' num2str(i)])) = subplot(3, 3, i);
    
    rdvec = 1+single(raw_data.head.idx.(matlab.lang.makeValidName(temp{i})));
    if nargin > 1
        rdvec = rdvec(subset_indices);
    end
    
    plot( rdvec ); title(temp{i}, 'Interpreter', 'none')
    hl = [hl, h.(matlab.lang.makeValidName(['x' num2str(i)]))];
end
linkaxes(hl, 'x');

end

function parpool_setup(requested_pp)
if nargin < 1
    % default pool num
    requested_pp = 16;
end

% check current instances
gcp_info = gcp('NoCreate');

% limit num workers to allowed number
pc_info = parcluster('local');
max_pp = pc_info.NumWorkers;
if requested_pp > max_pp; requested_pp = max_pp; end

% set-up parallel pool
if isempty(gcp_info.isvalid)
    
    parpool('local', requested_pp);
    
% else % <optional>
    %     % boost the number of workers if necessary
    %     numw = gcp_info.NumWorkers;
    %
    %     if numw < requested_pp
    %         delete(gcp('nocreate'));
    %         parpool('local', requested_pp);
    %     end
end
end

function [dmtx] = noise_adjust(nfile, d_MRD_h, data_dwell_time)
    % ====================
    % Load data 
    % ====================
    dsetin      = ismrmrd.Dataset(nfile, 'dataset');
    n_MRD_h     = ismrmrd.xml.deserialize(dsetin.readxml);
    noise_data  = dsetin.readAcquisition;
   
    disp(['Required Sens Map: ' d_MRD_h.measurementInformation.measurementDependency.measurementID ...
        ', Noise ID: ' n_MRD_h.measurementInformation.measurementID]);
    
    % ====================
    % Coil checks go here
    % ====================
    
    % ====================
    % Grab data
    % ====================
    
    n_samples = double(noise_data.head.number_of_samples(1));
    n_channels = double(noise_data.head.active_channels(1));
    
    % assuming Siemens using 2 averages:
    noise_ind = 256;
    nt2 = zeros(n_samples, noise_ind, n_channels);
    
    for i = 1:noise_ind
        nt2(:,i,:)=  reshape(noise_data.data{i}, [n_samples, 1, n_channels ]);
    end
    
    % Scale factor for receiver bandwidth (fixed for Siemens, empirically derived)
    rBW = n_MRD_h.acquisitionSystemInformation.relativeReceiverNoiseBandwidth;
        
    % ====================
    % Calc decorrelation matrix
    % ====================
    
    n_scaling = rBW * data_dwell_time / (noise_data.head.sample_time_us(1)*1e-6);
    dmtx = ismrm_calculate_noise_decorrelation_mtx(nt2, n_scaling ); 
    % figure,imagesc(abs(dmtx));

end
