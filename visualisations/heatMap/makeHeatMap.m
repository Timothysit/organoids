function makeHeatMap(spikeMatrix, option) 

    % outputs a heatmap figure handle for your spikes 
    % Assume numSamp x 60 matrix 
    
    % INTPUT
        % spikeMatrix 
            % number of samples x numChannels binary spike matrix 
            % where 0 means no spike and 1 means spike 
            % current implementation is based on 8 x 8 grid with 60
            % electrodes, and with each corner removed 
        % option 
            % default option is 'count', which shows the spike count 
            % if option is 'rate', then this outputs the firing rate
            % instead 
            % note the current implementation assumes a sampling rate of
            % 25kHz
         
    
    
    % Author: Tim Sit 
    % Last Update: 20180607
    
    if ~exist('option')
        option = 'count'; 
    end  
    
    
    spikeCount = sum(spikeMatrix); 
    
    if strcmp(option,'rate')
        spikeCount = spikeCount / (size(spikeMatrix, 1) / 25000);
    end 
    
    
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
    
    % flip things (organoid project)
    heatMatrix = heatMatrix'; 
    
    % make heatmap, whilst setting NA values to white
    h = imagesc(heatMatrix); 
    set(h, 'AlphaData', ~isnan(heatMatrix))
    
    
    % colormap(viridis)
    
    % attempt to make electrode with 0 spike (excluded) 
    % (not the same as NA, which are not recorded)
    % myColorMap = viridis(256); 
    % myColorMap(1, :) = 0; % 0 is black, 1 is white
    % colormap(myColorMap); 
   
    
    cb = colorbar;
    if strcmp(option, 'count') 
        ylabel(cb, 'Spike count')
    elseif strcmp(option, 'rate') 
        ylabel(cb, 'Spike rate (Hz)')
    end 
    cb.TickDirection = 'out';
    cb.Location = 'Southoutside';
    cb.Box = 'off';
    
    aesthetics
    removeAxis
    
    % make it square-ish
    set(gcf, 'Position', [100, 100, 800, 700])
    
    % set font size 
    set(gca, 'FontSize', 14)

end 