% Run through all .mat files within directory, extract spike matrix, 
% save in file with same name appendded with 'stats'

%% Options for spike detection 

method = 'Manuel'; % spike detection method 
multiplier = 5; % threshold multiplier (will depend on which method you chose) 

%% Loop through files in directory 

files = dir('*.mat'); 
progressbar
for file = 1:length(files)
    load(files(file).name, 'electrodeMatrix'); 
    data = electrodeMatrix;
    spikeTrain = getSpikeMatrix(data, method, multiplier);
    spikeTimes = findSpikes(spikeTrain);
    % save
    [filepath, name, ext] = fileparts(files(file).name);
    fileName = strcat(name, '_info', '.mat'); 
    save(fileName, 'spikeTimes', 'spikeTrain', '-v7.3');
    clear electrodeMatrix spikeTrain spikeTimes
    progressbar(file / length(files))
end