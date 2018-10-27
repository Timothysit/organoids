%% Batch network spike detection 
% sparse down sample function 
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/The_Mecp2_Project/feature_extraction/matlab/analysis_functions_ts/processSparse/'))

dirPath = '/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/spike_data/';
networkSpikeDetectMethod = 'simple'; 
networkSpikeDetectParam.minChannel = 7;
networkSpikeDetectParam.binWidth = 0.1; % in seconds
fs = 25000;
jitterCriteria = true; % remove simultaneous spikes

batchNetworkSpikeDetect(dirPath, networkSpikeDetectMethod, networkSpikeDetectParam, fs, jitterCriteria)

%% Batch network spike visualisation 
% aesthetics 
addpath(genpath('/media/timothysit/Seagate Expansion Drive1/TS_MEA_TOOLBOX/visualisation/'))

saveFolder = '/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/spike_data/networkSpikeData/networkSpikeData_minChannel7_binWidth100ms_jitterTrue/';
networkSpikeData = load('/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/spike_data/networkSpikeData/networkSpikeData_minChannel7_binWidth100ms_jitterTrue.mat');
networkSpikeData = networkSpikeData.networkSpikeData;
folderPath = '/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/data/all_mat_files';
batchVisualiseNetworkSpike(folderPath, networkSpikeData, 'multitrace', saveFolder)
fprintf('Finished making plots. \n')