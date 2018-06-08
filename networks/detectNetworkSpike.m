function [networkSpikeMatrix, networkSpikeTimes] = detectNetworkSpike(spikeMatrix, fs, networkSpikeDur, minChannel)
%DETECTNETWORKSPIKE Detect network spikes and gives back matrix summarising
% which electrodes or channels were participating in each of them 

% requires the downSampleSum function, which can be obtained here: 
% https://github.com/Timothysit/mecp2

% Author: Tim Sit 
% Context: For the Organoid Project
% Last Update: 20180518

% TODO: consider sliding window rather than binning

% INPUT 
    % spikeMatrix     | numSamp x numChannel spike matrix 
    % fs              | sampling frequency (Hz) 
    % networkSpikeDur | The duration of each network spike to search for
    % (seconds) 
% OUTPUT 
    % netSpikeMatrix  | numChannel x numNetworkSpike matrix 
    % where numNetworkSpike is the number of network spikes detected
    % where 1 indicate participation of a channel in that network spike

    % networkSpikeTimes | The time of the start of each networkSpike
    % (seconds) 

% STEP 1: Bin the spikeMatrix

newSampleFreq = networkSpikeDur * fs; % bin it to the duration of network spikes
recordingDuration = size(spikeMatrix, 1) / fs; 
newSampleNum = recordingDuration * newSampleFreq; 
% there should be a faster way to do this, but I can't think of it ATM

% TODO: return error if the above is not an integer, or forcibly round it
downSpikeMatrix = downSampleSum(spikeMatrix, newSampleNum);

% STEP 2: Convert any time with more than one spike into having just one (in order
% to use sum to identify the number of "active" electrodes)

downSpikeMatrix(downSpikeMatrix >= 1) = 1;

% STEP 3: Search for times where <minChannel> channels has more than one spike

networkSpikeFrames = find(sum(downSpikeMatrix, 2) >= minChannel); 
% sum along the second dimension, ie. the all columns (each column
% represent an electrode, so sum over all electrodes) 

% STEP 4: Generate networkSpikeMatrix and convert networSpikeFrames to
% seconds

networkSpikeMatrix = downSpikeMatrix(networkSpikeFrames, :); 

networkSpikeTimes = networkSpikeFrames / newSampleFreq;
    
    
end
