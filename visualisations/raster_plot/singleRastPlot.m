function singleRastPlot(spikeTrain) 
% single raster plot for one electrode / trace 
% assume spikeTrain to have dimensions numSamples x 1 
% where 1 in the spikeTrain represent yes spike, and 0 for no spike 

spikePos = find(spikeTrain == 1); 
plot([spikePos spikePos], [0 1], 'k')
ylim([-0.5 1.5]) 
xlim([1 length(spikeTrain)])
removeAxis()
aesthetics() 
lineThickness(1)
% TODO: add x-axis for time
end 
