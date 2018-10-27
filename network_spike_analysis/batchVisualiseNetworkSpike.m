function batchVisualiseNetworkSpike(folderPath, networkSpikeData, method, saveFolder)
%BATCHVISUALISENETWORKSPIKE Summary of this function goes here
%{
Goes through network spike data, plots and save network spikes. 

INPUT 
folderPath       | path to folder containing the raw .mat data (they are used to
plot the network spikes 

networkSpikeData | cell containing information about network spikes,
produced with batchNetworkSpikeDetect.m

method           | method to perform network spike visualisation 
                   'multitrace': put electrodes in rows (show only active
                   electrodes)
                   'grid'      : put electrodes in a grid (show all electrodes)

saveFolder       | Folder path to save the plotted figures



OUTPUT
saves image files of network spike to directory.

%}

%% Define file path and parameters

fileList = dir(strcat(folderPath, '*.mat'));
fileNames = {fileList(:).name}';
fileNames(ismember(fileNames,{'.','..'})) = [];% remove . and ..

fileWithNetworkSpike = find(cell2mat(networkSpikeData(:, 1)) > 0 );
numFileWithNetworkSpike = length(fileWithNetworkSpike);

binWidth = networkSpikeData{1, 4};
fs = 25000;

rawRootPath = '/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/data/all_mat_files/';


%% Define Filter 
% butterworth filter 
lowpass = 600; 
highpass = 8000; 
wn = [lowpass highpass] / (fs / 2); 
filterOrder = 3;
[b, a] = butter(filterOrder, wn); 
% filteredData = filtfilt(b, a, double(data)); 


%% Visualisation
vsep = 35; % vertical separation between traces
elecMax = 10; % maximum number of channels to plot
maxNetworkSpikeToPlot = 10; % maximum amount of network spikes to make for each file

for file = 1:numFileWithNetworkSpike 
    loadPath = strcat(rawRootPath, networkSpikeData{fileWithNetworkSpike(file), 7});
    rawData = load(strrep(loadPath{1}, '_thresholdSpikes', '')); %{1} basically does cell2str
    electrodeMatrix = rawData.dat;
    if size(electrodeMatrix, 1) < size(electrodeMatrix, 2)
        electrodeMatrix = electrodeMatrix';
    end 
    filteredMatrix = filtfilt(b, a, double(electrodeMatrix)); 
    partElectrodesList = networkSpikeData{fileWithNetworkSpike(file), 3};
    timeRangeList = networkSpikeData{fileWithNetworkSpike(file), 2};
    % cell containing the participating electrodes for that file 
    % will be further unpacked in the for loop to get vectors
    networkSpikeCount = 1;
    numNetworkSpike = networkSpikeData{fileWithNetworkSpike(file), 1};
    for networkSpike = 1:numNetworkSpike
        if strcmp(method, 'multitrace')
            try 
                networkSpikeStart = (timeRangeList(networkSpike) - binWidth/2) * fs;
                networkSpikeEnd = (timeRangeList(networkSpike) + binWidth/2) * fs;
            catch 
                fprintf('Network spike out of time range, skipping...\n')
                continue 
            end 
            timeRange = round(networkSpikeStart) : round(networkSpikeEnd);
            partElectrodes = partElectrodesList{networkSpike};
            multiElectrodeTrace(filteredMatrix, timeRange, partElectrodes, vsep, elecMax)
            saveName = strcat(saveFolder, networkSpikeData{fileWithNetworkSpike(file), 7});
            saveName = strrep(saveName, 'thresholdSpikes', strcat('spike_', num2str(networkSpikeCount)));
            saveName = strrep(saveName, '.mat', '.png'); 
            saveas(gcf, saveName{1});
            close all 
            networkSpikeCount = networkSpikeCount + 1;
        elseif strcmp(method, 'grid')
            % TODO: Add grid method for batch visualisation
        end
        if networkSpike > maxNetworkSpikeToPlot % exit early (prevent too many plots)
            fprintf('Too many network plots for this file, exiting early.\n')
            continue  
        end 
    end 
end 

end

