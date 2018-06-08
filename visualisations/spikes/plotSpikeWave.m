function  plotSpikeWave(trace, spikeTrain, method, fs, durationInSec)
%PLOTSPIKEWAVEFORM Plots the waveform of detected spikes
    % note that this assumes you already have the detected spikes 
    % this function does not detect spikes, use spikeDetect.m for that
%   INPUTS: 
    % trace      | raw or filtered trace from a single electrode. numSamp x 1
    % spikeTrain | binary vector containing 1 for spike locations, and 0 otherwise 
    % method     | method to plot detected waveforms, either just 'simple', or
    % aligned with the peaks 'peak' 
    % fs         | sampling frequency of your raw or filtered trace 
    % durationInSec | Duration (in seconds) in which you want to plot
    % with the spikes centred, usually 0.002 will be a good start. 
%   OUTPUTS: 
    % a figure of the waveform 

    [spikeWaves, averageSpikes] = spikeAlignment(trace, spikeTrain, fs, durationInSec); 
    plotSpikeAlignment(spikeWaves, method, fs, durationInSec);
% Author: Tim Sit 
% Last update: 20180427


    
end

