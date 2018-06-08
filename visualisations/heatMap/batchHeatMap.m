% loop through all the spike files in your directory, make a heatmap, then
% save

files = dir('*.mat'); 


for file = 1:length(files)
    data = load(files(file).name, 'cSpikes'); 
    spikes = data.cSpikes; 
    
    makeHeatMap(spikes);
    
    % set range (for MECP2 project only, rough guess is 0 - 100) 
    caxis(log10([10, 3500]))
    colorbar('southoutside')
    % save
    [filepath, name, ext] = fileparts(files(file).name);
    
    % set title of the heatmap 
    % title(name(1:end-5), 'Interpreter', 'none') % -5 to remove "_info"
    % disabled interpreter so that underscore doesn't become subscript
    
    % log scale tick marks
%     minSpike = log10(10); 
%     maxSpike = log10(3500); 
%     ticks = linspace(minSpike, maxSpike, 6);
%     YTickLabels = cellstr(num2str(round(ticks(:), 2), '10^%d'));
%     
%     cbh=colorbar('h');
%     set(cbh,'YTick',YTickLabels)
% this doens't work yet...


    % set position 
    set(gcf, 'Position', [100 100 900 800]);
    
    fileName = strcat(name(1:end-5), 'heatMap', '.png'); 
    saveas(gcf, fileName)
    clear spikes
    close all 
end