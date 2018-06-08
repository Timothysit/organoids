% prez's spike detection method 

figure 

%% Raw data 
subplot(100, 1, [1 40])
plot(data)
title('Raw data')
aesthetics()


%% elliptical filter, threshold and spikes
% mainly based on Wave_clus
% elliptical filter
% no idea what 0.1 and 40 do ...
par.detect_order = 4; % default, no idea why specifically this number
par.ref_ms = 1.5; % refractor period in ms, not sure when this is going to be used...
sr = 25000; % sampling rate 
fmin_detect = 300; 
fmax_detect = 8000;
[b,a] = ellip(par.detect_order,0.1,40,[fmin_detect fmax_detect]*2/sr);
% FiltFiltM does the same thing, but runs slightly faster
% the parameters are default found in wave_clus
filteredData = filtfilt(b, a, double(data)); 

% finding threshold and spikes
med = median(filteredData); 
s = std(filteredData); 
multiplier = 5;
threshold = med - multiplier*s; 
spikeTrain = filteredData < threshold; % negative threshold
spikeTrain = double(spikeTrain);
numSpike = sum(spikeTrain);

% raster plot 
subplot(100, 1, [48 51])
singleRastPlot(spikeTrain) 

subplot(100, 1, [56 100])
plot(filteredData); 
title('4th order Elliptical Filter 300 - 8000Hz, 5SD threshold')
% threshold 
hold on; 
plot([1 length(data)], [threshold threshold], '-')
aesthetics()

%% some diagonistics 

% number of psikes 
% calculate an  average firing rate from that






