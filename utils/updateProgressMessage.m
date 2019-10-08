function updateProgressMessage(n, ntot, tbase, len, freq)

%% Defaults
if nargin<5 || isempty(freq)
    freq = 1; % freq of full message updates; else just print '.' on each call
end

if nargin<4 || isempty(len)
    len = 80; % padded string length
end

if nargin<3 || isempty(tbase)
    t = toc;
else
    t = toc(tbase);
end


%% Make string

if mod(n, freq)==1
    % update times
    secPerN = t/n;
    tRemEstimate = secPerN * (ntot-n)/60;
    % clear previous message
    if n>1
        fprintf(repmat('\b',1,len + freq-1));
    end
    % update message
    msg = sprintf('\nfinished %i of %i.  (%2.2f min elapsed, ~%2.2f min remaining; %2.2f s/ea)',n, ntot, t/60, tRemEstimate, secPerN);
    fprintf(pad(msg, len, '-'));
else
    fprintf('.')
end
    
end


