function [spikeMatrix, filteredMatrix] = getSpikeMatrix(data, method, multiplier, L)
% Loop through detectspikes for each channel 
% Assume input matrix: numSamp x numChannels

% L is the loss ratio parameter, for 'cwt' method only 

if ~exist('L', 'var')
    L = 0; % no bias; cost of spike = cost of noise
end

spikeMatrix = zeros(size(data)); 
filteredMatrix = zeros(size(data)); 
for j = 1:size(data, 2)
    [spikeMatrix(:, j), finalData, ~] = detectSpikes(data(:, j), method, multiplier, L);
    % fprintf(num2str(j)) % this was for debugging, to see which electrode 
    % is making an error
    filteredMatrix(:, j) = finalData;
end 

end 