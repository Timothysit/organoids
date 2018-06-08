function adjM = getAdjM(spikes, method, downSample, lag)
% getAdjM computes the weighted adjacency matrix of given spike counts 

% INPUT 
    % spikes: spike counts, expect dimensions to be numSample x numChannel
    % method: string argument for the method used to compute adjacency 
        % correlation (OKAY)
        % partial correlation (WORKING)
        % xcorr : correlation with lag (TODO)
        % mutual information (TODO)
    % downsample: the new total number of samples to downsample to
    % set to 0 if you don't want it to be downsampled
    % this may be required for partial correlation, or xcov, as it may
    % exceed matlab array size limit 
    % lag: the time lag to do xcorr
    
if ~exist('lag', 'var') 
    lag = 1; % should be 15ms as per Schoeter et al 2015
end 

% unsparse the matrix if it is sparse

% author: Tim Sit 
% Last Update: 20180302
% Pending tasks: 
% 1) coherence 
% 2) mutual information, if that is of interest 

% some implementation notes specific to the mecp2 project 
% 12 minute recording, 25kHz 
% may need to downsample to about 250Hz for xcov to work (array size of
% about 8gb)

% for partial correlation, the recommended sampling frequency is: 


if issparse(spikes) 
    spikes = full(spikes); 
end 

% downsample

if downSample ~= 0
   spikes = downSampleSum(spikes, downSample);
end 


% calculate weighted adjacency matrix 
        
if strcmp(method, 'correlation') 
    adjM = corr(spikes); 
elseif strcmp(method, 'partialcorr')
    adjM = partialcorr(spikes);
    % [Pxcorr] = Plagged(spikes); 
elseif strcmp(method, 'xcorr')
    spikes = full(spikes);
    c = xcorr(spikes,lag);
    c(lag + 1, :) = [ ]; % remove middle; remove 0 lag 
    cMax = max(c); % take max across all range of time lag 
    numChannel = size(spikes, 2); % reshape this so it resembles an adjacency matrix 
    adjM = reshape(cMax, numChannel, numChannel); 
elseif strcmp(method, 'mutInfo')
    % mutual information 
    
elseif strcmp(method, 'tileCoef')
    fs = 25000;
    numChannel = size(spikes, 2);
    dtv = lag;
    Time = [0, length(spikes) / fs]; 
    combChannel = nchoosek(1:numChannel, 2); 
    adjM = NaN(numChannel, numChannel); 
    for channelPairNum = 1:size(combChannel, 1)
        electrode_1 = spikes(:, combChannel(channelPairNum, 1)); 
        electrode_2 = spikes(:, combChannel(channelPairNum, 2));
        N1v = sum(electrode_1); 
        N2v = sum(electrode_2); 
        spike_times_1 = find(electrode_1 >= 1) / fs; 
        spike_times_2 = find(electrode_2 >= 1) / fs;
        tilingCoef = sttc(N1v, N2v, dtv, Time, spike_times_1, spike_times_2);
        adjM(combChannel(channelPairNum, 1), combChannel(channelPairNum, 2)) = tilingCoef; 
        adjM(combChannel(channelPairNum, 2), combChannel(channelPairNum, 1)) = tilingCoef; 
        % this computes only the triangle, so I have to also copy that to the
        % other side        
    end
    % assign diagonal values to 1 
    % adjM(logical(eye(size(adjM)))) = 1;
    
    % only assign diagonal value of 1 to active electrodes, otherwise, NaN
    activeElectrode = find(sum(spikes >= 1)); 
    for atvE = 1:length(activeElectrode) 
        adjM(activeElectrode(atvE), activeElectrode(atvE)) = 1; 
    end 
    
end 


%% Analysis on raw / filtered signal, ie. coherence etc. 






end 