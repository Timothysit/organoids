function coefVal = spikePartCoef(networkSpikePair)
% spikePartCoef computes the network spike participation coefficient
% between a pair of channels/nodes/electrodes. This coefficient summarises
% the extent to which two neurons/nodes/electrodes participate in the same
% network spike

% Author: Tim Sit 
% Last Update: 20180519 

% note that the input networkSpikePair can be obtained from the
% networkSpikeMatrix, which can be geneated with the detectNetworkSpike.m
% function (which in turn depends on having the spikeMatrix)

% INPUT 
    % networkSpikePair  | 2 x N matrix, where N is the total number of
    % network spikes observed in the NETWORK (not just the pair). Each row 
    % represent one channel 
    
% OUTPUT 
    % coefVal           | The network spike participation coeficient 
    % value interpretation
    % -1 : fully anticorrelated network spike participation
    % 0  : there is a balance between anticorrelated participation and correlated participation 
    % 1  : fully correlated network spike participiation 

%% Check that there is at least one network spike in each channel 

if sum(networkSpikePair(1, :)) < 1 || sum(networkSpikePair(2, :)) < 1
    coefVal = NaN; 
    return 
end 

%% Find N, the network spikes where either channel participate in 
% Step 1: Find the total number of non-overlapping network spikes 
% participated by the node $x_i$ and $x_j$, this will be $N$.
% In other words, the number of total network spike
% participated by either node will be called $N$. 

participation = sum(networkSpikePair, 1); 
N_spikePair = networkSpikePair(:, participation >= 1); 
N = size(N_spikePair, 2);

%% Compute the participation scores and get the coefficient
% Step 2: For each network spike $n$ in $N$, compute the "participation score" according to the table below

%| channel 1| channel 2 | participation score |
% |-------|-------|---------------------|
% |     0 |     1 |                  -1 |
% |     1 |     0 |                  -1 |
% |     1 |     1 |                   1 |

participation_score = sum(N_spikePair, 1); 
% if the sum is 1, then it must be either {0, 1} or {1, 0} 
% if the sum is 2, then it must be {1, 1} 

participation_score(participation_score == 1) = -1; 
participation_score(participation_score == 2) = 1;

% Step 3: Sum the scores over all $N$ network spikes, 
% and divide that by the value $N$ 
coefVal = sum(participation_score) / N; 

end 