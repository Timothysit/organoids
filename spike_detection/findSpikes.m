function spikeTimes = findSpikes(spikeTrain) 
% goes through a spike train / spike matrix and find spikes 
% assume 1 = spike and no spikes otherwise 
% if spike matrix, assume dimension to be samples x electrodes 
% output is a cell array of vectors containing the spike times in frames

% according to this, a for loop works about as well as a vectorised
% solution
%https://www.mathworks.com/matlabcentral/answers/229711-using-find-on-each-column-of-a-matrix-independently

spikeTimes = cell(1, size(spikeTrain, 2)); % pre-allocate

for n = 1:size(spikeTrain, 2)
    spikeTimes{n} = find(spikeTrain(:, n) == 1); 
end

end 