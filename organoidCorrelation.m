%% Organoid Project - Correlational Analysis 

% author: Tim Sit 
% Last Update: 20180508
% For The Organoid Project; data recorded by Dr. Susanna Mierau and Stefano
% assume data loaded and spikes extraced using the main script 
% organoidProject.m

%% Calculation of correlation
% choice of time bin will be critical here!
spikes = a;
method = 'correlation'; 
downSample = 0; 
lag = 1;
adjM = getAdjM(spikes, method, downSample, lag);

%% Spike tiling coefficient method to find correlation 

% find spike times (in seconds) 

electrode_1 = spikeMatrix(:, 19); 
electrode_2 = spikeMatrix(:, 20);

N1v = sum(electrode_1);  
N2v = sum(electrode_2);
Time = [1 / fs, length(spikeMatrix) / fs]; 
spike_times_1 = find(electrode_1 == 1) / fs; 
spike_times_2 = find(electrode_2 == 1) / fs;

% test out for a pair of spike trains
dtv = 0.01; % in seconds, the delay time to look for conincidental spikes 
tileCoef = sttc(N1v, N2v, dtv, Time, spike_times_1, spike_times_2); 

%% tiling coefficient over all electrode pairs 
% this will take a while...
method = 'tileCoef';
downSample = 0;
lag = 0.04; 
adjM = getAdjM(spikeMatrix, method, downSample, lag); 



%% Plot Adacenecy matrix 
figure
h = imagesc(adjM);

% make NaN values white 
% set(h, 'AlphaData', ~isnan(adjM))
yticks([1, 10:10:60])
xticks([1, 10:10:60])
aesthetics 
set(gca,'TickDir','out'); 
cb = colorbar;
ylabel(cb, 'Correlation')
cb.TickDirection = 'out';
caxis([0 1])
cb.Ticks = 0:0.1:1; 
set(gca,'TickDir','out'); 
cb.Location = 'Southoutside';
cb.Box = 'off';
set(gca, 'FontSize', 14)

yLength = 800; 
xLength = yLength * 1.2; 
set(gcf, 'Position', [100 100 xLength yLength])


%% Network plot! 

% need to emphasise on making this asetheticly pleasing
figure
goodElectrodes = 1:60; 
% TODO: make sure electrode 31 has no correlation, 
% TODO: make sure locations are correct
plotAdj(adjM, goodElectrodes')

yLength = 800; 
xLength = yLength * 1.2; 
set(gcf, 'Position', [100 100 xLength yLength])


%% Relationship between correlation coefficient and distance 

% set up distance matrix 
gridLength = 8;
electrodePairs = npermutek(1:gridLength, 2);

% remove the 4 corners 
corners = [1; 8; 57; 64]; 
electrodePairs(corners, :) = [ ];


% distanceMatrix = pdist(electrodePairs, 'euclidean'); % nope
distanceMatrix = distmat(electrodePairs);

% remove self connections 
distanceMatrix(logical(eye(size(distanceMatrix)))) = NaN;
corrMatrix = adjM;
corrMatrix(logical(eye(size(corrMatrix)))) = NaN;

% plot relationship between distance and corr coef 
figure 
dotSize = 60;
elecSpacingCoef = 200;
scatter(distanceMatrix(:) * elecSpacingCoef, corrMatrix(:), dotSize, 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0); 
xlabel('Distance (\mum)') 
ylabel('Tiling coefficient') 
set(gca, 'FontSize', 14)
set(gca,'TickDir','out'); 

% plot individual electrodes with different colours 
figure 
dotSize = 60;
elecSpacingCoef = 200; % they are 200 micro-meters apart (\mu m)
for elec = 1:length(channels) 
    scatter(distanceMatrix(:, elec) * elecSpacingCoef, corrMatrix(:, elec), dotSize, 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0); 
    hold on 
end 
xlabel('Distance (\mum)')
ylabel('Tiling coefficient') 
set(gca, 'FontSize', 14)
set(gca,'TickDir','out');

%% Network spike participation 

networkSpikeDur = 10 * 10^-3; % 10 ms
minChannel = 5; 
fs = 25000;
[networkSpikeMatrix, networkSpikeTimes] = detectNetworkSpike(spikeMatrix, fs, networkSpikeDur, minChannel);

% visualise participation 
figure
imagesc(networkSpikeMatrix')
colormap(flipud(gray))  % black and white, with 1 being black
aesthetics
ylim([1 size(spikeMatrix, 2)])
ax = gca;
ax.YTick = [1 5:5:60];
ylabel('Electrode number')
xlabel('Network spike number') 
set(gca, 'FontSize', 14)
set(gca,'TickDir','out')

%% Network spike participation correlation 

% requres the networkSpikeMatrix 
networkSpikeDur = 5 * 10^-3; % 10 ms
minChannel = 5; 
fs = 25000;
[networkSpikeMatrix, networkSpikeTimes] = detectNetworkSpike(spikeMatrix, fs, networkSpikeDur, minChannel);

numChannel = size(spikeMatrix, 2);
combChannel = nchoosek(1:numChannel, 2); 
adjM = NaN(numChannel, numChannel); 
    for channelPairNum = 1:size(combChannel, 1)
        electrode_1 = networkSpikeMatrix(:, combChannel(channelPairNum, 1))'; 
        electrode_2 = networkSpikeMatrix(:, combChannel(channelPairNum, 2))';
        coefVal = spikePartCoef([electrode_1; electrode_2]); 
        adjM(combChannel(channelPairNum, 1), combChannel(channelPairNum, 2)) = coefVal; 
        adjM(combChannel(channelPairNum, 2), combChannel(channelPairNum, 1)) = coefVal; 
    end
% assign diagonal values to 1, but only if the participated in at least one
% network spike 
% adjM(logical(eye(size(adjM)))) = 1;

activeElectrode = find(sum(networkSpikeMatrix >= 1)); 
% activeElectrode = find(sum(spikeMatrix >= 1)); 
for atvE = 1:length(activeElectrode) 
    adjM(activeElectrode(atvE), activeElectrode(atvE)) = 1; 
end 

%% Visualise network spike participation for each electrode 

figure
imagesc(networkSpikeMatrix')
colormap(flipud(gray))  % black and white, with 1 being black
aesthetics
ylim([1 size(spikeMatrix, 2)])
ax = gca;
ax.YTick = [1 5:5:60];
ylabel('Electrode number')
xlabel('Network spike number') 
set(gca, 'FontSize', 14)
set(gca,'TickDir','out')

%% Visualise participation correlation 

% adjacency matrix 
figure
h = imagesc(adjM);

% make NaN values white 
set(h, 'AlphaData', ~isnan(adjM))


aesthetics 
set(gca,'TickDir','out'); 
cb = colorbar;
ylabel(cb, 'Correlation')
cb.TickDirection = 'out';
% cb.Ticks = 0:0.1:1; 
set(gca,'TickDir','out'); 
cb.Location = 'Southoutside';
cb.Box = 'off';
set(gca, 'FontSize', 14)


% network plot

figure
goodElectrodes = 1:60; 
% TODO: make sure electrode 31 has no correlation, 
% TODO: make sure locations are correct
plotAdj(adjM, goodElectrodes')


