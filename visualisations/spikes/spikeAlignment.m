function [spikeWaves, averageSpike] = spikeAlignment(data, spikeTrain, fs, durationInSec)
% aligns the spikes by giving you the voltage signal just before, during,
% and after the spike. This allows you to visualise the detected spike to
% verify them. 

% Written by: Tim Sit 2017-2018
% Last Update: 20180427 

% INPUT 
    % data : the voltage-trace / any form of continuous time series 
    % expect it to be a vector of n x 1, where n is the number of samples 
    
    % spikeTrain : a binary vector of size also n x 1 
    % 0 means there are no spike in that sample, 1 means spike 
    
    % fs : sampling frequency in Hz 
    
    % durationInSec : the duration of the time series to be extracted for each
    % spike. The spike will be in the centre of this extracted window.
    % Input in seconds. 
    % for microelectrode arrays, 2 * 10^-3 seconds is usually sufficient
    
% OUTPUT 
    % spikeWaves : a numS x duration matrix
    % where numS is the number of spikes detected 
    % and each row represent the time series, with the spike being detected
    % in the middle of that time sreies 
    
% TODO 
    % set default variable values

% fs = 25000; 
% duration = fs * 1 * 10^-3; % assume each spike is one millisec 
% duration = fs * 2 * 10^-3;
duration = fs * durationInSec; 

% what we need is to find spike from the spike train 
% if it is a spike, then we plot duration frames around that point 
% then we hold, and search for the next spike

% subplot(1, 2, 1)
spikeTimes = find(spikeTrain == 1);
spikeStore = zeros(length(spikeTimes), 2*round(duration/2)+1); % pre-allocation is paramount here
for i = 1:length(spikeTimes)
     spikePoint = spikeTimes(i);
     % spikeRaw = data(spikePoint - round(duration/2):(spikePoint + round(duration/2))); 
     spikeStore(i, :) = data(spikePoint - round(duration/2):(spikePoint + round(duration/2))); 
     % spikeStore(i, :) = data(spikePoint - floor(duration/2):(spikePoint + floor(duration/2)));
end
% aesthetics()

averageSpike = mean(spikeStore);
spikeWaves = spikeStore; 
%subplot(1, 2, 2)
%title('Average spike')
%plot(averageSpike)
%aesthetics()

end