

function chanMaps = loadChanMaps()


ksroot = fileparts(fileparts(mfilename('fullpath')));
chanMapFiles = dir(fullfile(ksroot, 'configFiles', '*.mat'));

idx = 1;
chanMaps = [];
orig_state = warning;
warning('off','all')

for c = 1:numel(chanMapFiles)
    if strcmp(chanMapFiles(c).name(1),'.')
        % skip hidden metadata files
        continue
    else
        q = load(fullfile(ksroot, 'configFiles', chanMapFiles(c).name));
        
        cm = createValidChanMap(q, chanMapFiles(c).name);
        if ~isempty(cm)
            if idx==1; chanMaps = cm; else; chanMaps(idx) = cm; end
            idx = idx+1;
        end
    end
end

warning(orig_state);
