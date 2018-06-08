%% The organoid Project 
% Author: Tim Sit 
% Last Update: 20180426

%% Load data 

% spontaneous activity / media / TTX % 20180413
cd('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/data/mat_files/')

% ChR2: 20180503
cd('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/data/organoid_20180503_light/mat_files')

% 20180508 files - ChR2
cd('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/data/organoid_20180508/mat_files')

% 20180518 files - Friday 
cd('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/data/organoid_20180518/mat_files')


%% Load some MEA data processing code 

% TODO: rebrand those code into some sort of MEA processing code
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/analysis_functions_ts/'))
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/organoid_data_analysis/'))

% scale bar 
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/chenxinfeng4-scalebar-4ca920b/'))
% heatmaps 
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/heatMap/'))
% human colours 
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/XKCD_RGB/'))
% cwt spike detection 
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/continuous_wavlet_transform/'))

% orgaoind project 
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project')); 

%% Set up some parameters

fs = 25000; % sampling frequency

%% Filter data 
% note that this is already done within spike detection functions


%%  grid trace 
electrodeMatrix = dat'; 
figure 
gridTraceVtwo(electrodeMatrix, 1000) % note that it is rotated for mecp2/organoid project  (compared to HPC datasets)

%% specific electrode plot 
% this can actually be done with a polished makeTrace function
electrodeNum = 6;
singleTrace = filteredMatrix(:, electrodeNum);

figure
plot(singleTrace);
xlim([1 length(electrodeMatrix)])
aesthetics

yLength = 350; 
xLength = yLength * 5;
set(gcf, 'Position', [100 100 xLength yLength])

% adjust y and x axis to make a bit of room for my scale bar 
yOffset = 1;
ylim([min(electrodeMatrix(:, electrodeNum)) - yOffset, max(electrodeMatrix(:, electrodeNum))])

xOffset = 100000; 
xlim([1 - xOffset, length(singleTrace)])

removeAxis 
sb = scalebar; 
sb.YLen = 5; 
sb.XLen = 1000000; 
sb.YUnit = '\muV';
sb.XUnit = 'seconds';
sb.Position = [-100000 -744.5];

%% Convert X-axis to seconds to comapre with raster plot 
timeBins = 5; % 5 second separation between marks
timePoints = 1:floor(length(electrodeMatrix) / (fs * timeBins)); 
xticks(timePoints * fs * timeBins); 
xticklabels(string(timePoints * timeBins));

%% spike detection (entire grid)
method = 'Manuel';
% method = 'cwt';
multiplier = 6;
L = 0; 
timeRange = 1: 600 * fs;
[spikeMatrix, filteredMatrix] = getSpikeMatrix(electrodeMatrix(timeRange, :), method, multiplier, L);

%% Positive threhsold for artefact removal 

% spikeMatrix(filteredMatrix < -100) = 0; % doesn't really work

artefactLoc = find(filteredMatrix < - 100);
removeDur = 0.2; % seconds to remove around the centre of the artefact

for aLoc = 1:length(artefactLoc) 
    spikeMatrix(artefactLoc(aLoc) - removeDur * fs : artefactLoc(aLoc) + removeDur * fs) = 0;
end


%% visualise number of spikes in a heatmap 
figure
makeHeatMap(spikeMatrix)
% note that makeHeatMap is modified for organoid project (tranposed) 

%% Truncate some of the spike matrix  

% only take the first 175 seconds for slice 3 
% 20180413 slice 5: 225 seconds
% 20180413 slice 4_baseline: 280 seconds
% 20180503 slice 2 recording 3: 360 seconds
% 20180503 slice 3 recording 3: 360 seconds
% 20180503 slice 3 recording 4: 90 seconds
% 20180503 slice 4 recording 1: 90 seconds
% 20180503 slice 6 recording 1: 620 seconds
timeRange = 1: 360 * fs; % in seconds
spikeMatrix_cut = spikeMatrix(timeRange, :); 
filteredMatrix_cut = filteredMatrix(timeRange, :);


% spikeMatrix = spikeMatrix(1:225 *25000, :);

%% plot spikes over time for each electrode

% TODO: make function to do this 
% down sample spike matrix to sum them in time bins 
recordDuration = length(spikeMatrix) / fs;
downSpikeMatrix = downSampleSum(spikeMatrix, recordDuration * 1/5); % 5 second time bins
figure
plot(downSpikeMatrix)
aesthetics
set(gca,'TickDir','out'); 
lineThickness(2) 
ylabel('Spike count') 
xlabel('Time bin (5 seconds)')
set(gca, 'FontSize', 14)


%% Raster plot




%% Heatmap raster plot! 

figure
recordDuration = length(spikeMatrix) / fs;
downSpikeMatrix = downSampleSum(spikeMatrix, recordDuration * 1/5); 

% Delete certan time points and replace with NA 
% For media drop / TTX drop purpose 
% deleteTime = [21 22, 38,39]; % note that these are in time bins for slice
% 4 
% deleteTime = [23 24]; % TTX drop time for slice 5
% downSpikeMatrix(deleteTime, :) = NaN; 

% h = imagesc(downSpikeMatrix'); 

% alternative raster plot that uses frequency rather than spike count

h = imagesc(downSpikeMatrix' ./5); 

aesthetics 
ylabel('Electrode') 
xlabel('Time (s)')
cb = colorbar;
% ylabel(cb, 'Spike count')
ylabel(cb, 'Spike rate (Hz)') 
cb.TickDirection = 'out';
% cb.Ticks = 0:5; % for slice 5 specifically
set(gca,'TickDir','out'); 
cb.Location = 'Southoutside';
cb.Box = 'off';
set(gca, 'FontSize', 14)


set(h, 'AlphaData', ~isnan(downSpikeMatrix')) % for NaN values

timeBins = 5; % 5 second separation between marks
% timePoints = 1:timeBins:floor(length(spikeMatrix) / fs); 
timePoints = 0:5:floor(length(spikeMatrix) / fs); 
yticks([1, 10:10:60])
xticks(timePoints); 
xticklabels(string(timePoints * 5));
% xticklabels(string(timePoints -1 ));

yLength = 800; 
xLength = yLength * 2; 
set(gcf, 'Position', [100 100 xLength yLength])



%% Within Burst Heatmap raster plot 

figure
recordDuration = length(spikeMatrix) / fs;
newFs = 100; % new frequency to sample, 
newSampFreq = recordDuration * newFs; % turn it into a 100Hz 
downSpikeMatrix = downSampleSum(spikeMatrix, newSampFreq);

startTime = 236; % in seconds
endTime = 238; % in seconds
h = imagesc(downSpikeMatrix(startTime * newFs : endTime *newFs, :)'); 
% h = imagesc(downSpikeMatrix')
aesthetics 
ylabel('Electrode') 
% xlabel('Time (s)')
cb = colorbar;
ylabel(cb, 'Spike count')
cb.TickDirection = 'out';
% cb.Ticks = 0:1; % for slice 5 specifically
set(gca,'TickDir','out'); 
cb.Location = 'Southoutside';
cb.Box = 'off';
set(gca, 'FontSize', 14)

% set(h, 'AlphaData', ~isnan(downSpikeMatrix')) % for NaN values

% timeBins = 5; % 5 second separation between marks
timePoints = 0:50:550; 
yticks([1, 10:10:60])
xticks(timePoints); 
xticklabels(string(timePoints / 1000));

yLength = 800; 
xLength = yLength * 2; 
set(gcf, 'Position', [100 100 xLength yLength])



%% plot detected spike waveforms 
figure
timeRange = 135 * fs : 185 * fs; % use 1:length(spikeMatrix) for all
% timeRange = 1:length(spikeMatrix); 
electrodeNum = 45;
trace = filteredMatrix(timeRange, electrodeNum);
spikeTrain = spikeMatrix(timeRange, electrodeNum); 
durationInSec = 0.009; 
% plotSpikeWave(trace, spikeTrain, 'peak', fs, durationInSec)
[spikeWaves, averageSpikes] = spikeAlignment(trace, spikeTrain, fs, durationInSec); 
plotSpikeAlignment(spikeWaves, 'peak', fs, 0.008); % note that input here has to be shorter than durationInSec
aesthetics 
lineThickness(2)
removeAxis 
yLength = 300; 
xLength = yLength * 21/9;
set(gcf, 'Position', [100 100 xLength yLength])
sb = scalebar; 
sb.YLen = 20;
sb.YUnit = 'a.u.';
sb.XUnit = 'ms';
sb.XLen = 20;
sb.Position = [2 -60];
xlim([1, 200])

%% Version 2 of plot detected spike waveform 
% the visual effect I want to go for here is transparent individual spikes 
% and then a thick mean waveform
figure
% timeRange = 135 * fs : 185 * fs; % use 1:length(spikeMatrix) for all
timeRange = 1:length(spikeMatrix); 
electrodeNum = 48;
trace = filteredMatrix(timeRange, electrodeNum);
spikeTrain = spikeMatrix(timeRange, electrodeNum); 
durationInSec = 0.009; 
% plotSpikeWave(trace, spikeTrain, 'peak', fs, durationInSec)
[spikeWaves, averageSpikes] = spikeAlignment(trace, spikeTrain, fs, durationInSec); 
plotSpikeAlignment(spikeWaves, 'peakghost', fs, 0.008); % note that input here has to be shorter than durationInSec
aesthetics 
lineThickness(2)
removeAxis 
yLength = 300; 
xLength = yLength * 21/9;
set(gcf, 'Position', [100 100 xLength yLength])
sb = scalebar; 
sb.YLen = 10;
sb.YUnit = '\muV';
sb.XUnit = 'ms';
sb.XLen = 20; % note that this is in frames, need to manually convert to ms
sb.Position = [2 -20];
xlim([1, 200])
ylim([-35, 20])

% direct export eps to make the transparency work 
print(gcf, '-opengl','-depsc', '-r600', '/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/avespike_20180518_slice7_record1_electrode34.eps')

% note that if you use the menu to export as eps, transparency won't work.
% (ie. the mean trace won't be placed on top of the gray traces properly
% for some reason)



%% Superimpose electrodes 
startTime = 105; 
endTime = 110; 
figure
for electrode = 1:length(channels) 
    timeRange = startTime * fs : endTime * fs; % use 1:end for all
    trace = filteredMatrix(timeRange, electrode);
    spikeTrain = spikeMatrix(timeRange, electrode); 
    durationInSec = 0.02; 
    % plotSpikeWave(trace, spikeTrain, 'peak', fs, durationInSec)
    [spikeWaves, averageSpikes] = spikeAlignment(trace, spikeTrain, fs, durationInSec); 
    plotSpikeAlignment(spikeWaves, 'peak', fs, 0.008);
    hold on 
end 
aesthetics 
lineThickness(2)
removeAxis 


%% Plot detected spike waveforms in 3D
figure
plotSpikeAlignment(subset_spikeWaves, 'peak3D', fs, 0.008); 
aesthetics 
lineThickness(2)
removeAxis

%% Remove some spikes, then plot again 

% a simple threshold 
threshold = 100;
subset_spikeWaves = spikeWaves(~any(spikeWaves > threshold, 2), :); % search rows 

%% Tune threshold until no spikes in no MEA case 

% an alternative is just get a ball-park value, eg. testing increments of
% 0.5 on the multiplier 

%% STATS: Spike after media drop, then TTX (for slice 4 only)

% slice 4 

% plot the number of spikes under different time windows 
% for each electrode 

% number of spikes before media for ea. electrode 

startSpikes = sum(spikeMatrix(1:100*fs, :));
startRate = startSpikes / 100; 

% number of spikes after media 

mediaSpikes = sum(spikeMatrix(110*fs : 185 *fs, :)); 
mediaRate = mediaSpikes / 75;

% number of spikes after TTX

ttxSpikes = sum(spikeMatrix(195 * fs : 300 *fs, :)); 
ttxRate = ttxSpikes / 105;

sliceFour(:,1) = startRate; 
sliceFour(:, 2) = mediaRate; 
sliceFour(:, 3) = ttxRate;

% get it into a format that works more nicely with R 
electrodeCount = 1;
for i = 1:3:length(sliceFour) * 3 
    sliceFourDf(i).electrodeName = electrodeCount;
    sliceFourDf(i).condition = 'Control';
    sliceFourDf(i).spikeRate = sliceFour(electrodeCount, 1); 
    sliceFourDf(i+1).electrodeName = electrodeCount;
    sliceFourDf(i+1).condition = 'Media';
    sliceFourDf(i+1).spikeRate = sliceFour(electrodeCount, 2); 
    sliceFourDf(i+2).electrodeName = electrodeCount;
    sliceFourDf(i+2).condition = 'TTX';
    sliceFourDf(i+2).spikeRate = sliceFour(electrodeCount, 3); 
    electrodeCount = electrodeCount + 1;
end 

save('sliceFourDf', 'sliceFourDf')

%% Plot light times on top of the raster plot 
figure
subplot(10, 1, 1:2) 

light = zeros(720 * fs, 1); 
% specify light on times
light(120 * fs : 180 *fs) = 1; 
light(240 * fs : 300 * fs) = 1; 
light(360 * fs: 420 * fs) = 1;
light(480 * fs : 540 * fs) = 1; 
light(600 * fs : 660 * fs) = 1; 

plot(light)
aesthetics
removeAxis 
lineThickness(2)


subplot(10, 1, 3:10) 
downSpikeVec = downSampleSum(spikeMatrix(:, 49), 720 / 5);
bb = bar(downSpikeVec);
bb.FaceColor = 'black';
bb.EdgeAlpha = 0; 
xlim([1, length(downSpikeVec)])
removeAxis

%% Plot raster plot on top of the filtered trace 
figure
subplot(10, 1, 1:2) 
% h = imagesc(downSpikeMatrix(:, 22)');  % heatmap approach 
% timeRange = 110.5 * fs: 111.5 * fs;
% spikeVec = fullSpikeMatrix(timeRange, 22); % slice 5 

spikeVec = spikeMatrix(:, 45); 
% singleRastPlot(spikeVec)
% [spikeVec, filteredVec] = detectSpikes(electrodeMatrix(1:225*fs, 33), method, multiplier, L);
% singleRastPlot(spikeVec)

recordDuration = length(spikeVec) / fs;
downSpikeVec = downSampleSum(spikeVec, recordDuration * 1/2);
h = imagesc(downSpikeVec'); 
removeAxis
cb = colorbar;
cb.TickDirection = 'out';
ylabel(cb, 'Spike count')
cb.Box = 'off';
set(gca, 'FontSize', 14)

% bar chart approach 
% bb = bar(downSpikeVec);
% bb.FaceColor = 'black';
% bb.EdgeAlpha = 0; 
% xlim([1, length(downSpikeVec)])
% removeAxis

subplot(10, 1, 3:10) 
% filteredVec = filteredMatrix(timeRange, 22);
filteredVec = filteredMatrix(:, 45);
plot(filteredVec); % slice 5 
xlim([1, length(filteredVec)]) % have no idea why matlab doesn't do this by default
aesthetics
removeAxis 
% scalebar


%% Combine pre and post TTX rasterplots 

% preTTXspike : slice 5 1 -115 sec 
% postTTXspike : slice 5 145 - 225 sec
% TTXNaN : the middle period (115 - 145) are all NaN values 

TTXNaN = NaN(30 *fs, 60); 

fullSpikeMatrix = [preTTXspike; TTXNaN; postTTXspike];


%% Multiple electrode trace plot within a burst
% specify time range of burst 
timeRange = 236 * fs : 238 *fs; 
% find the electrodes participating in a burst 
partElectrodes = find(sum(spikeMatrix(timeRange, :)) > 0);

% subplot method
figure 
for elec = 1:10 
subplot(10, 1, elec) 
plot(filteredMatrix(timeRange, partElectrodes(elec)))
aesthetics 
removeAxis
end 

% try offset method 
figure 
offSet = 0; 
for elec = 1:10
    plot(filteredMatrix(timeRange, partElectrodes(elec)) - offSet)
    hold on 
    offSet = offSet - 35;
end 
aesthetics 
removeAxis 

%% Single trace paper specific parameters 

% 20180518 slice 7 recording 1 
% electrode 34 
figure
plot(filteredMatrix(:, 34))
aesthetics
removeAxis
yLength = 500; 
xLength = yLength * 21/9;
set(gcf, 'Position', [100 100 xLength yLength])
ylim([-25, 25]) 
xlim([482 * fs, 487 * fs])
sb = scalebar;
sb.YUnit = '\muV'; 
sb.XUnit = 'ms'; 
sb.Position = [12060000 -20];

% electrode 41
figure
plot(filteredMatrix(:, 41))
xlim([275 *fs, 280 *fs])
ylim([-25, 25]) 
aesthetics 
removeAxis 
sb = scalebar;
sb.YUnit = '\muV'; 
sb.XUnit = 'ms'; 
yLength = 500; 
xLength = yLength * 21/9;
set(gcf, 'Position', [100 100 xLength yLength])
sb.Position = [6875000 -20];


% electrode 43 
figure 
plot(filteredMatrix(:, 43)) 
xlim([410 * fs, 415 * fs])
ylim([-25, 25]) 
aesthetics 
removeAxis 
sb = scalebar;
sb.YUnit = '\muV'; 
sb.XUnit = 'ms'; 
yLength = 500; 
xLength = yLength * 21/9;
set(gcf, 'Position', [100 100 xLength yLength])
sb.Position = [10252000 -20];

% electrode 48 

figure 
plot(filteredMatrix(:, 48)) 
xlim([318 * fs, 323 * fs])
ylim([-25, 25]) 
aesthetics 
removeAxis 
sb = scalebar;
sb.YUnit = '\muV'; 
sb.XUnit = 'ms'; 
yLength = 500; 
xLength = yLength * 21/9;
set(gcf, 'Position', [100 100 xLength yLength])
sb.Position = [318 *fs -20];

%% For paper: network burst grid plot 

burstStart = 110.75; 
burstEnd = 111.75;
burstMatrix = filteredMatrix(burstStart * fs : burstEnd * fs, :); 

figure 
gridTrace(burstMatrix, 1)
linkaxes

% let's try tight subplot for the gridTrace

figure 
gridTrace(burstMatrix, 1, [], 'tight', 1)
linkaxes

yLength = 600; 
xLength = yLength * 1.5;
set(gcf, 'Position', [100 100 xLength yLength])

% add a scalebar 
sb = scalebar; 
sb.YUnit = '\muV'; 
sb.XUnit = 'ms'; 
sb.YLen = 40; 
sb.XLen = 4000;



% plot individual traces 
figure 
plot(filteredMatrix(burstStart * fs : burstEnd * fs, 2))
xlim([1, length(burstMatrix)])
ylim([-100, 150])
aesthetics 
removeAxis 
yLength = 400; 
xLength = yLength * 21 / 9;
set(gcf, 'Position', [100 100 xLength yLength])
sb = scalebar; 
sb.YUnit = '\muV'; 
sb.XUnit = 'ms'; 
sb.YLen = 40; 
sb.Position = [100, -75];


%% Paper figure: raster plot 

% 0413 slice 6 
% note that electrode 24 needs to be 'grounded', there is abrupt change
% in baseline level 
% also note that the artefact around the 80 second mark needs to be removed
% artefact across all electrodes beteween 76.95 sec and 76.96 sec
% therefore remove spikes from all electrodes at that time

% the timing is quite important since there are spike like events happening
% within the 5 second as well... conincidence??? I think I will have to
% include them in principle, but I am skeptical.

spikeMatrix(:, 24) = 0;

removeSpikeElectrode = 1:60; 
removeSpikeElectrode([8, 47, 48, 49]) = [];
artefactStart = round(76.7 * fs); 
artefactEnd = round(77.2 * fs); 

spikeMatrix(artefactStart:artefactEnd, :) = 0;

%% Paper figure: pre and post-TTX drop 

% data: 0413 slice 4 drop then TTX 
% media drop time: 100 - 110 seconds
% TTX drop time: 185 - 190 seconds (but can make it 185 - 195 seconds
% Therefore, I think it will be good to make it 50 pre-post drop (it better
% shows the spontaneous activity before the pre-, if I select something
% like 30, there will be less spontaneous activity shown 

% Also have to make the time of media drop 0 (I think I will just use
% xticklabels to do that)


% remove TTX drop period: 

spikeMatrix(185 * fs: 195 * fs, :)  = NaN; 

% truncate the matrix to 135 sec to 245 sec 

spikeMatrix = spikeMatrix(135 * fs +1 : 245 * fs, :); 

% specific tick marks for the raster plot
timePoints = 1:2:150; 
xticks(timePoints); 
% xticklabels(string(timePoints / 1000));

xlab = {'-45', '-35', '-25', '-15', '-5', '', '5', '15', '25', '35', '45'};
xticklabels(xlab)
% specific dimensions for the raster plot 
yLength = 800; 
xLength = yLength * 1.2; 
set(gcf, 'Position', [100 100 xLength yLength])

% heatmap before and after 

preTTXspikeMatrix = spikeMatrix(1: 50 * fs, :);
postTTXspikeMatrix = spikeMatrix(60 * fs + 1: 110 * fs, :);

% turn nan values to zeros 
preTTXspikeMatrix(isnan(preTTXspikeMatrix)) = 0;


% pre-TTX heatmap
figure
makeHeatMap(preTTXspikeMatrix, 'rate')
set(gcf, 'Position', [100, 100, 800, 800 * 1])

% post-TTX heatmap

figure 
makeHeatMap(postTTXspikeMatrix, 'rate')
caxis([0, 0.3])
set(gcf, 'Position', [100, 100, 800, 800 * 1])

%% Paper figures: bar charts (using gramm)
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/piermorel-gramm-682ec28/'));
% bar chart of active electrodes before and after

activeElectrodePreTTX = sum(sum(preTTXspikeMatrix) > 0);
activeElectrodePostTTX = sum(sum(postTTXspikeMatrix) > 0);

g = gramm('x', {'pre-TTX', 'post-TTX'}, 'y', [activeElectrodePreTTX, activeElectrodePostTTX]);
g.set_names('x', '', 'y', 'Number of active electrodes')
g.geom_bar();
figure('Position',[100 100 800 800]);
g.set_text_options('base_size', 14)
g.axe_property('TickDir','out')
g.set_order_options('x', 0) % change the order of the bar plot 
g.set_color_options('map','d3_10')
g.draw(); 
set(gcf, 'Position', [100, 100, 300, 300 * 16/9])

% bar chart of firing frequency before and after 

ratePreTTX = sum(preTTXspikeMatrix) / (length(preTTXspikeMatrix) / fs); 
ratePreTTX = ratePreTTX(ratePreTTX > 0);
ratePostTTX = zeros(size(ratePreTTX)); 

sem = [std(ratePreTTX) / sqrt(length(ratePreTTX)), 0];

g = gramm('x', {'pre-TTX', 'post-TTX'}, 'y', [mean(ratePreTTX), mean(ratePostTTX)], ... 
    'ymax', [mean(ratePreTTX), mean(ratePostTTX)] + sem, ... 
    'ymin', [mean(ratePreTTX), mean(ratePostTTX)] - sem);
g.set_names('x', '', 'y', 'Spike rate (Hz)')
g.geom_bar();
figure('Position',[100 100 800 800]);
g.set_text_options('base_size', 14)
g.axe_property('TickDir','out', 'YLim', [0 0.1])
g.set_order_options('x', 0) % change the order of the bar plot 
g.set_color_options('map','d3_10')

% error bar 
g.geom_interval('geom','black_errorbar','dodge',0.8,'width',1);

g.draw(); 
set(gcf, 'Position', [100, 100, 300, 300 * 16/9])



%% Paper figure: spike detection threshold

% this is repurposed from compareSpikeDetect.m in the mecp2 project
electrodeNum = 45; 
timeRange = 1: 185 * fs - 1; % use 1:length(spikeMatrix) for all
data = electrodeMatrix(timeRange, electrodeNum);

[spikeTrain, finalData, threshold] = detectSpikes(data, 'Manuel', 6);
figure
ax1 = subplot(100, 1, [1 20]);
singleRastPlot(spikeTrain) 
% numSpike = sum(spikeTrain);
% title(['Manuel: Butterworth, mean - 5SD, 2.0ms RP,' s num2str(numSpike) s 'spikes'])
ax2 = subplot(100, 1, [21 100]); 
plot(finalData); 
% threshold 
hold on; 
plot([1 length(data)], [threshold threshold], '-')
aesthetics()
removeAxis()
xlim([1 length(data)])
linkaxes([ax1, ax2], 'x')

% TIME WINDOW TO FOCUS ON 
windowRange = [120 *fs, 170 *fs];
xlim(windowRange)

set(gcf, 'Position', [100, 100, 400 * 21/9, 400])

sb = scalebar;
sb.YUnit = '\muV'; 
sb.XUnit = 'ms'; 
sb.YLen = 20; 
sb.Position = [100 * fs, -75]; 

%% Paper figure: single spike resolution detection 
% note the key here is to adjust the xlimit 
% NOT the time range to do spike deteciotn, because then you may mess up
% the mean voltage level

% this is repurposed from compareSpikeDetect.m in the mecp2 project
% 0413 slice 4
electrodeNum = 45; 
timeRange = round(1) : round(300 * fs); % use 1:length(spikeMatrix) for all
data = electrodeMatrix(timeRange, 45);

[spikeTrain, finalData, threshold] = detectSpikes(data, 'Manuel', 5);
figure
% ax1 = subplot(100, 1, [1 20]);
% singleRastPlot(spikeTrain) 
numSpike = sum(spikeTrain);
% title(['Manuel: Butterworth, mean - 5SD, 2.0ms RP,' s num2str(numSpike) s 'spikes'])
% ax2 = subplot(100, 1, [21 100]); 
plot(finalData); 
% threshold 
hold on; 
plot([1 length(data)], [threshold threshold], '-')
aesthetics()
removeAxis()
xlim([135.09 * fs, 135.104 * fs])
% linkaxes([ax1, ax2], 'x')

set(gcf, 'Position', [100, 100, 400 * 21/9, 400])

sb = scalebar;
sb.YUnit = '\muV'; 
sb.XUnit = 'ms'; 
sb.YLen = 20; 
sb.XLen = 25;
sb.Position = [135.09*fs, -55]; 



