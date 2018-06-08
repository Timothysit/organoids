function adjM = dynamicAdjacency(spikes, numChunk, method)
%DYNAMICADJACENCY Computes the time varying adjacency matrix

%%% Input: 
      % spikes : numSamp x numChannel spike matrix. Assume sparse matrix
           % 1 = spike, 0 = no spike. 
           % TODO: implement full() so matrix type won't matter
      % numChunk : how many chunks to divide the matrix
      % method : method for calculating the adjacency matrix 
%%% Output
    % adjM : numChannl x numChannel matrix, not weighted.
    % wM   : weighted matrix (TODO)

windowLength = size(spikes, 1) / numChunk;  
chunks = [0:numChunk] .* windowLength;
chunks(1) = 1; 

adjM = zeros(size(spikes, 2), size(spikes, 2), numChunk); 

for n = 1:numChunk
    
    spikeChunk = spikes(chunks(n):chunks(n+1), :);  
    
    if strcmp(method, 'covariance') 
        covM = sparseCov(spikeChunk); 
        adjM(:, :, n) = covM > 1.25 * mean(covM); % thresholds
    end 

end 

adjM = double(adjM);

end

