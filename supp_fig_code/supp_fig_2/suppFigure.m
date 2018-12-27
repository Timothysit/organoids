%{
Code to generate the latest version of the supp. figures. 
If unspecified here, then code in organoidProject.m are the most updated
versions
Last update: 2018127
%}

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

% figure2eps
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/figure2epsV1-3'));


%% Update to supp figure 2 k and l 
% 20181227: Add dots on the bar plot (geom_jitter)
% data: 0413 slice 4 drop then TTX 
% media drop time: 100 - 110 seconds
% TTX drop time: 185 - 190 seconds (but can make it 185 - 195 seconds
% Therefore, I think it will be good to make it 50 pre-post drop (it better
% shows the spontaneous activity before the pre-, if I select something
% like 30, there will be less spontaneous activity shown 

% Also have to make the time of media drop 0 (I think I will just use
% xticklabels to do that)

% load data 
load('data/organoid_20180413/Organoid_180413_slice_4_MEA_3D_uncoated_drop_then_TTX.mat')
 
% spike extraction 
electrodeMatrix = dat'; 
method = 'Manuel'; 
multiplier = 6; % usual value: 5.5
L = 0; 
[spikeMatrix, filteredMatrix] = getSpikeMatrix(electrodeMatrix, method, multiplier, L);

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

% Paper figures: bar charts (using gramm)
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
g.set_names('x', '', 'y', 'Spike Frequency (Hz)')
g.geom_bar();
figure('Position',[100 100 800 800]);
g.set_text_options('base_size', 14)
g.axe_property('TickDir','out', 'YLim', [-0.005, 0.1])
g.set_order_options('x', 0) % change the order of the bar plot 
g.set_color_options('map','d3_10')

% error bar 
g.geom_interval('geom','black_errorbar','dodge',0.8,'width',1);

g.draw(); 


set(gcf, 'Position', [100, 100, 300, 300 * 16/9])


%% Attempt to make geom_jitter plot and a box plot 
condition = cell(1, length(ratePreTTX) * 2);
condition(1:length(ratePreTTX)) = {'pre-TTX'}; 
condition(length(ratePreTTX)+1:end) = {'post-TTX'}; 
gScatter = gramm('x', condition, 'y', [ratePreTTX, ratePostTTX]);
gScatter.geom_jitter();
gScatter.geom_bar();
gScatter.draw()

