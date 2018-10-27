function batchNetworkSpikeDetect(folderPath, method, networkSpikeDetectParam, fs, jitterCriteria)
%{ 
batchNetworkSpikeDetect performs network spike detection in batch 
and saves a table summarising the results

% Author: Tim Sit 
% Last update: 20181027 

% LOG 
% 20181027: Adding jitter method (in progress)

INPUT 
folderPath     | Path containing the spike files (.mat) you want to analyse 
method         | Method for performing network spike detection
jitterCriteria | If true, remove simultaneous spikes (potential artefact)
before looking for network spikes
fs             | sampling frequency


SPIKE DETECTION METHODS (<variable_to_provide_in_networkSPikeDetectParam>)

'binCount': simply bin the spike matrix of a given length <binWidth> (in
seconds) then look for time bin where the number of active channels 
(with at least one spike) is more than <minChannel>



OUTPUT

resultTable | table summarising the number of network spikes in each file,
and contains a cell containing specific information about those spikes


%}

fileList = dir(strcat(folderPath, '*.mat'));
fileNames = {fileList(:).name}';
fileNames(ismember(fileNames,{'.','..'})) = [];% remove . and ..

networkSpikeData = cell(length(fileNames), 6);

for file = 1:length(fileList) 
    data = load(strcat(folderPath, fileNames{file}));
    spikeMatrix = data.m5p5Spikes; % numSample x numChannel sparse double
    
    if jitterCriteria == true
        % jitter method adds a criteria whereby network spikes must be
        % non-simultaneous
        % an improvement to this is to add a parameter jitterLength controlling how
        % much (in samples) the spikes must be apart.
        eaFrameSpikeCount = sum(spikeMatrix, 2); % sum all channel per frame
        simultaneousSpikeIndex = find(eaFrameSpikeCount > 1);
        % turn simultaneous spike frames to have no spikes
        spikeMatrix(simultaneousSpikeIndex, :) = 0;
    end 
    
    if strcmp(method, 'simple')
        newSampleRate = 1 / networkSpikeDetectParam.binWidth;
        newSampleNum = size(spikeMatrix, 1) / fs * newSampleRate;
        downSpikeMatrix = sparseDownSample(spikeMatrix, newSampleNum, 'sum'); 
        % downSampleSum or sparseDownSample
        
        % binarise each time bin so that it is only spike/no spike
        spikeBinary = downSpikeMatrix > 0;
        channelSumBinary = sum(spikeBinary, 2);
        networkSpikeFrames = find(channelSumBinary > networkSpikeDetectParam.minChannel);
        networkSpikeTimes = networkSpikeFrames / newSampleRate;
        numNetworkSpike = length(networkSpikeTimes);
        
        networkSpikeMatrix = downSpikeMatrix(networkSpikeFrames, :);
        participatingElectrodes = cell(size(networkSpikeMatrix, 1), 1);
        
        for nSpike = 1:size(networkSpikeMatrix, 1)
            participatingElectrodes{nSpike} = find(networkSpikeMatrix(nSpike, :) > 0);
        end 
        
        % participatingElectrodes = find(channelSumBinary == 1); % still in
        % development, need to convert to matrix actually...
        
        networkSpikeData{file, 1} = numNetworkSpike;
        networkSpikeData{file, 2} = networkSpikeTimes;
        networkSpikeData{file, 3} = participatingElectrodes;
        networkSpikeData{file, 4} = networkSpikeDetectParam.binWidth;
        networkSpikeData{file, 5} = networkSpikeDetectParam.minChannel;
        networkSpikeData{file, 6} = 'simple';
        networkSpikeData{file, 7} = fileNames(file);        
    end 
end 

save('networkSpikeData', 'networkSpikeData')
end 