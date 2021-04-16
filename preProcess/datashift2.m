function rez = datashift2(rez, do_correction)

rez.iorig = 1:rez.temp.Nbatch;

if  getOr(rez.ops, 'nblocks', 1)==0
    return;
end

% display status of flag for whether to apply upsampled datashift (==1), datashift rounded to interger channel spacing (==2),
% or no correction at all (just makes the datashift plots in that case)
do_correction

% flag for additional datashift plot of corrected data
debugPlot =  getOr(ops, 'fig', 1)==2;

ops = rez.ops;

% The min and max of the y and x ranges of the channels
ymin = min(rez.yc);
ymax = max(rez.yc);
xmin = min(rez.xc);
xmax = max(rez.xc);

% Determine the average vertical spacing between channels. 
% Usually all the vertical spacings are the same, i.e. on Neuropixels probes. 
dmin = median(diff(unique(rez.yc)));
fprintf('Computing datashift alignment...\n')
fprintf('\tvertical pitch size is %d \n', dmin)
rez.ops.dmin = dmin;
rez.ops.yup = ymin:dmin/2:ymax; % centers of the upsampled y positions

% Determine the template spacings along the x dimension
% dminx = median(diff(unique(rez.xc)));
yunq = unique(rez.yc);
mxc = zeros(numel(yunq), 1);
for j = 1:numel(yunq)
    xc = rez.xc(rez.yc==yunq(j));
    if numel(xc)>1
       mxc(j) = median(diff(sort(xc))); 
    end
end
dminx = median(mxc);
fprintf('\thorizontal pitch size is %d \n', dminx)

rez.ops.dminx = dminx;
nx = round((xmax-xmin) / (dminx/2)) + 1;
rez.ops.xup = linspace(xmin, xmax, nx); % centers of the upsampled x positions


if  getOr(rez.ops, 'nblocks', 1)==0
    rez.iorig = 1:rez.temp.Nbatch;
    return;
end



% binning width across Y (um)
dd = 10;% 5;
%   - rest of drift computation relies on indices of this sampling res.
%   - If dd res is too low/high for electrode spacing, drift fit will be suboptimal
% min and max for the range of depths
dmin = ymin - 1;
dmax  = 1 + ceil((ymax-dmin)/dd);


spkTh = 10; %ops.ThPre; %10;% [10] same as the usual "template amplitude", but for the generic templates

% Extract all the spikes across the recording that are captured by the
% generic templates. Very few real spikes are missed in this way. 
[st3, rez] = standalone_detector(rez, spkTh);
%%

% detected depths
dep = st3(:,2);

% min and max for the range of depths
dmin = ymin - 1;
dep = dep - dmin;

dmax  = 1 + ceil(max(dep)/dd);
Nbatches      = rez.temp.Nbatch;
batchSec      = ops.NT/ops.fs;

% which batch each spike is coming from
batch_id = st3(:,5); %ceil(st3(:,1)/dt);

% preallocate matrix of counts with 20 bins, spaced logarithmically
F = zeros(dmax, 20, Nbatches);
for t = 1:Nbatches
    % find spikes in this batch
    ix = find(batch_id==t);
    
    % subtract offset
    dep = st3(ix,2) - dmin;
    
    % amplitude bin relative to the minimum possible value
    amp = log10(min(99, st3(ix,3))) - log10(spkTh);
    
    % normalization by maximum possible value
    amp = amp / (log10(100) - log10(spkTh));
    
    % multiply by 20 to distribute a [0,1] variable into 20 bins
    % sparse is very useful here to do this binning quickly
    M = sparse(ceil(dep/dd), ceil(1e-5 + amp * 20), ones(numel(ix), 1), dmax, 20);    
    
    % the counts themselves are taken on a logarithmic scale (some neurons
    % fire too much!)
    F(:, :, t) = log2(1+M);
end

%%
% determine registration offsets
ysamp = dmin + dd * [1:dmax] - dd/2;
[imin,yblk, F0, F0m] = align_block2(F, ysamp, ops.nblocks);

if isfield(rez, 'F0')
    d0 = align_pairs(rez.F0, F0);
    % concatenate the shifts
    imin = imin - d0;
end

%%
if getOr(ops, 'fig', 1)  
    figure;
    set(gcf, 'Color', 'w')
    
    % plot the shift trace in um
    plot(imin * dd)
    box off
    xlabel('batch number')
    ylabel('drift (um)')
    title('Estimated drift traces')
    drawnow
    
    % raster plot of all spikes at their original depths
    ii = st3(:,3)>=spkTh;
    st_depth0 = st3(ii,2);
    st_depthD = st_depth0 + imin(batch_id(:)) * dd;
    xs = (.5:1:Nbatches)*batchSec;
    
    H = figure;
    set(H, 'name', 'driftMap')
    hax1 = subplot(2+debugPlot,1,1); hold on; box off
    hax2 = subplot(2+debugPlot,1,2); hold on; box off

    for j = spkTh:100
        % for each amplitude bin, plot all the spikes of that size in the
        % same shade of gray
        ix = st3(:, 3)==j; % the amplitudes are rounded to integers
        thisCol = [1 1 1] * max(0, 1-j/60);
        plot(hax1, st3(ix, 1)/ops.fs, st_depth0(ix), '.', 'color', thisCol) % the marker color here has been carefully tuned
        plot(hax2, st3(ix, 1)/ops.fs, st_depthD(ix), '.', 'color', thisCol) % the marker color here has been carefully tuned
    end
    axis tight
    linkaxes([hax1,hax2]);
    plot(hax1, xs, imin * dd +diff(ylim(hax2))/2, '-r');

    title(hax1, 'Drift map: Initial', 'interp','none')
    ylabel(hax1, 'Init spike position (um)')
    title(hax2, 'Drift map: Corrected', 'interp','none')    
    ylabel(hax2, 'Shifted spike position (um)')
    xlabel(hax2, 'time (sec)')

    % name figure
    addFigInfo(ops, H);

end
%%
% convert to um 
dshift = imin * dd;

% NO@!!@  Please stop messing with time for drift correction!!!
% % this is not really used any more, should get taken out eventually
% [~, rez.iorig] = sort(mean(dshift, 2));

if do_correction
    if do_correction==2
        % round shifts to integers of channel spacing
        dmin = median(diff(unique(rez.yc)));
        dshift = round(dshift./dmin).*dmin;
    end
    % sigma for the Gaussian process smoothing
    sig = rez.ops.sig;
    % register the data batch by batch
    dprev = gpuArray.zeros(ops.ntbuff,ops.Nchan, 'single');
    for ibatch = 1:Nbatches
        dprev = shift_batch_on_disk2(rez, ibatch, dshift(ibatch, :), yblk, sig, dprev);
    end
    fprintf('time %2.2f, Shifted up/down %d batches. \n', toc, Nbatches)
else
    fprintf('time %2.2f, Skipped shifting %d batches. \n', toc, Nbatches)
end
% keep track of dshift 
rez.dshift = dshift;
% keep track of original spikes
rez.st0 = st3;

rez.F = F;
rez.F0 = F0;
rez.F0m = F0m;

% next, we can just run a normal spike sorter, like Kilosort1, and forget about the transformation that has happened in here 

% additional debug plotting
% - slow & just for confirmation of what the ACTUAL saved whitened data looks like
% - read back the newly shifted data and plot driftmap
% - should look just like shifted driftmap, but with clipping on any channels that were shifted out of dat file register

if debugPlot && any(rez.dshift)
    
    %% One more extraction & plot to show actual shifts applied
    % This is slow, but reassuring.  ....disable when certain of parameters
    
    spkTh = 10;%ops.ThPre; %10;% [10] same as the usual "template amplitude", but for the generic templates
    
    % Extract all the spikes across the recording that are captured by the
    % generic templates. Very few real spikes are missed in this way.
    % - don't pass out rez struct this time
    % - we don't want/need any of this to carry over, just checking our work
    [st3] = standalone_detector(rez, spkTh);
    
    % detected depths
    dep = st3(:,2);
    
    % min and max for the range of depths
    dmin = ymin - 1;
    dep = dep - dmin;
    
    dmax  = 1 + ceil(max(dep)/dd);
    Nbatches      = rez.temp.Nbatch;
    batchSec      = ops.NT/ops.fs;
    
    % which batch each spike is coming from
    batch_id = st3(:,5); %ceil(st3(:,1)/dt);
    
    % preallocate matrix of counts with 20 bins, spaced logarithmically
    F = zeros(dmax, 20, Nbatches);
    for t = 1:Nbatches
        % find spikes in this batch
        ix = find(batch_id==t);
        
        % subtract offset
        dep = st3(ix,2) - dmin;
        
        % amplitude bin relative to the minimum possible value
        amp = log10(min(99, st3(ix,3))) - log10(spkTh);
        
        % normalization by maximum possible value
        amp = amp / (log10(100) - log10(spkTh));
        
        % multiply by 20 to distribute a [0,1] variable into 20 bins
        % sparse is very useful here to do this binning quickly
        M = sparse(ceil(dep/dd), ceil(1e-5 + amp * 20), ones(numel(ix), 1), dmax, 20);
        
        % the counts themselves are taken on a logarithmic scale (some neurons
        % fire too much!)
        F(:, :, t) = log2(1+M);
    end
    
    figure(H)
    % raster plot of all spikes at their original depths
    ii = st3(:,3)>=spkTh;
    st_depth0 = st3(ii,2);
    
    hax3 = subplot(3,1,3); hold on; box off

    for j = spkTh:100
        % for each amplitude bin, plot all the spikes of that size in the
        % same shade of gray
        ix = st3(:, 3)==j; % the amplitudes are rounded to integers
        thisCol = [1 1 1] * max(0, 1-j/60);
        plot(hax3, st3(ix, 1)/ops.fs, st_depth0(ix), '.', 'color', thisCol) % the marker color here has been carefully tuned
    end
    axis tight
    linkaxes([hax1,hax2,hax3]);
    plot(hax3, xs, rez.dshift +diff(ylim(hax2))/2, '-r');

    title(hax3, 'Drift map: Applied', 'interp','none')    
    ylabel(hax3, 'Drifted spike position (um)')
    xlabel(hax3, 'time (sec)')
end


end % main function

