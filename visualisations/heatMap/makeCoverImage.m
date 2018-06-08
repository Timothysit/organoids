% make cover image for my thesis 

% it should contain the spike heat map, then the correlation heatmap, then
% the control heatmap 

%% load data 



%% do the calculations 
% choose something DIV28 so that there is full heatmap image
% tc043 DIV28 might be somewhere to start 
load('/media/timothysit/Seagate Expansion Drive/The_Mecp2_Project/feature_extraction/matlab/data/schroter_spikes/HP_tc043_DIV28_spikes.mat')

spikes = m5Spikes;

% some pre-processing 
samplingRate = 25000;
recordDuration = 9 * 60; % in seconds
spikes = spikes(end - 60 * samplingRate * 9 + 1:end, :); % get last 9 minutes  
downFreq = 25;
downSpikes = sparseDownSample(spikes, recordDuration * downFreq, 'sum');
noSpikeElectrode = find(sum(downSpikes) == 0); 
downSpikes(:, noSpikeElectrode) = []; 

% minSpikes = 10; 
% maxSpikes = 10^6; % maximum threshold just in case something exceed biological range
% % TODO: Logical indexing instead of find
% goodElectrodes = find(sum(spikes) >= minSpikes); 
% spikes = spikes(:, goodElectrodes); 
% numActiveElectrode = length(goodElectrodes); % for the record

% metrics 
adjM = getAdjM(downSpikes, 'correlation', 0, 0); 
aveControl = ave_control(adjM); 

%% coercse the control heatmap to go properly 
controlMetric = aveControl;

% set control metric to percentage of total control 
controlMetric = controlMetric ./ sum(controlMetric) * 100;

%%%%%%%%%%%%%%%%%%%%%% Make heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
heatMatrix = zeros(8, 8); 
heatMatrix(1, 1) = NaN; 
heatMatrix(1, 8) = NaN; 
heatMatrix(8, 1) = NaN;
heatMatrix(8, 8) = NaN; 

% this assumes all electrode present
if length(controlMetric) == 60
    heatMatrix(2:7) = controlMetric(1:6);
    heatMatrix(9:56) = controlMetric(7:54); 
    heatMatrix(58:63) = controlMetric(55:60);
elseif length(controlMetric) == 59
    % basically, grounded electrode 5
    heatMatrix(5) = NaN; 
    heatMatrix(2:4) = controlMetric(1:3); 
    heatMatrix(6:7) = controlMetric(4:5); 
    heatMatrix(9:56) = controlMetric(6:53);
    heatMatrix(58:63) = controlMetric(54:59);
elseif length(controlMetric) == 58
    % basically, grounded electrode 5 and 16
    heatMatrix(5) = NaN;
    heatMatrix(16) = NaN;
    heatMatrix(2:4) = controlMetric(1:3); 
    heatMatrix(6:7) = controlMetric(4:5);
    heatMatrix(9:15) = controlMetric(6:12); 
    heatMatrix(17:56) = controlMetric(13:52); 
    heatMatrix(58:63) = controlMetric(53:58);
end 




%% make figure 

figure 
subplot(1, 3, 1) 
makeHeatMap(spikes) 

subplot(1, 3, 2) 
imagesc(adjM)
removeAxis

subplot(1, 3, 3) 
% make heatmap, whilst setting NA values to white
h = imagesc(heatMatrix); 
set(h, 'AlphaData', ~isnan(heatMatrix))
removeAxis


