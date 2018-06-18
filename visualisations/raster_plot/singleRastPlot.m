function singleRastPlot(spikeTrain, option) 
% single raster plot for one electrode / trace 
% assume spikeTrain to have dimensions numSamples x 1 
% where 1 in the spikeTrain represent yes spike, and 0 for no spike 

% option : 
    % 'line' raster plot with straight line at each spike time
    % 'dot' raster plot with a circular dot at each spike time 

if ~exist('option')
        option = 'line'; 
end 
    

spikePos = find(spikeTrain == 1); 
if strcmp(option, 'line')
    plot([spikePos spikePos], [0 1], 'k')
elseif strcmp(option, 'dot')
    plot([spikePos spikePos], [0 0], 'k', 'Marker', '.', 'MarkerSize', 10)
end


ylim([-0.5 1.5])
xlim([1 length(spikeTrain)])
removeAxis()
aesthetics() 
lineThickness(1)
% TODO: add x-axis for time
end 
