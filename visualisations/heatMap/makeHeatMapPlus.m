function makeHeatMapPlus(spikeMatrix) 
    % INCOMPLETE, last update: 20180127
    % same as makeHeatMap, but with added electrode trace on the side
    % outputs a heatmap figure handle for your spikes 
    % Assume 60 x 60 matrix 
    spikeCount = sum(spikeMatrix); 
    
    % log(spikeCount) % log scale 
    
    % reshape it to display properly 
    % note that this part is quite hard-coded, don't use for general
    % heatmap, only use the imagesc part for non-MEA data
    
    heatMatrix = zeros(8, 8); 
    heatMatrix(1, 1) = NaN; 
    heatMatrix(1, 8) = NaN; 
    heatMatrix(8, 1) = NaN;
    heatMatrix(8, 8) = NaN; 
    
    % this assumes all electrode present
    if length(spikeCount) == 60
        heatMatrix(2:7) = spikeCount(1:6);
        heatMatrix(9:56) = spikeCount(7:54); 
        heatMatrix(58:63) = spikeCount(55:60);
    elseif length(spikeCount) == 59
        % basically, grounded electrode 5
        heatMatrix(5) = NaN; 
        heatMatrix(2:4) = spikeCount(1:3); 
        heatMatrix(6:7) = spikeCount(4:5); 
        heatMatrix(9:56) = spikeCount(6:53);
        heatMatrix(58:63) = spikeCount(54:59);
    elseif length(spikeCount) == 58
        % basically, grounded electrode 5 and 16
        heatMatrix(5) = NaN;
        heatMatrix(16) = NaN;
        heatMatrix(2:4) = spikeCount(1:3); 
        heatMatrix(6:7) = spikeCount(4:5);
        heatMatrix(9:15) = spikeCount(6:12); 
        heatMatrix(17:56) = spikeCount(13:52); 
        heatMatrix(58:63) = spikeCount(53:58);
    end 
    
    % not sure what to do if its not 60 x 60...
    
    % make heatmap, whilst setting NA values to white
    subplot(1, 2, 1)
    h = imagesc(heatMatrix); 
    set(h, 'AlphaData', ~isnan(heatMatrix))    
    colormap(viridis)
    cb = colorbar;
    ylabel(cb, 'Spike count')
    
    subplot(1, 2, 2)
    gridTrace(spikeMatrix, 1000)
    
    aesthetics
    removeAxis
    
    % make it square-ish
    set(gcf, 'Position', [100, 100, 800, 700])

end 