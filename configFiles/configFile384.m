ops.chanMap             = '';
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if no chanMap file

% sample rate
ops.fs = 40000;  

% frequency for high pass filtering (150)
ops.fshigh = 300;%150;   

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0;%0.1; 

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 4];  

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = 10;  

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.95; 

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/100; % [def = 1/50]  % unk. implementation of this param, but firing rate cutoffs should up to user post-hoc

% number of samples to average over (annealed from first to second value) 
ops.momentum = [20 400]; 

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 70; % [def = 30] % must be wider for U-probe spacing

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8; 

%% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction. 

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
% ops.Nfilt               = 1024; % max number of clusters
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
% Define batches in seconds of data, not abstract bit chunks
    batchSec                = 5;  % define batch number in seconds of data [def = 2] ...longer if <100 channels, else batch reordering gets screwy
    batchSkips             = ceil(60/batchSec); % do high-level assessments at least once every second of data
    ops.NT                  = ceil(batchSec*ops.fs/32)*32 +ops.ntbuff; % convert to 32 count increments of samples
% ops.NT                  = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = batchSkips; % compute whitening matrix from every N-th batch
ops.nskip               = batchSkips;  % how many batches to skip for determining spike PCs
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available

%%