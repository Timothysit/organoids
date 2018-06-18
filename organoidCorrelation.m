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

% x and y axis labels 
xlabel('Electrode') 
ylabel('Electrode')

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
xLength = yLength * 0.90; 
set(gcf, 'Position', [100 100 xLength yLength])


%% Network plot! 

% need to emphasise on making this asetheticly pleasing
figure
goodElectrodes = 1:60; 
% TODO: make sure electrode 31 has no correlation, 
% TODO: make sure locations are correct
plotAdj(adjM, goodElectrodes')

yLength = 800; 
xLength = yLength * 1.1; 
set(gcf, 'Position', [100 100 xLength yLength])

% remember to export it ussing export rather than save as
% or else the aspect ratio will not be correct

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
% ylabel('Tiling coeffecient')
ylabel('Correlation')
set(gca, 'FontSize', 14)
set(gca,'TickDir','out'); 
aesthetics
ylim([0 1])
set(gcf, 'Position', [100, 100, 500 * 16/9, 500])

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

%% Paper figure: distance distribution plot 

% aim of this is to show that we get highly correlated nodes across the
% grid (near and far)

% one way to do this it to plot a histogram of the distances for
% connections with "high correlation", eg. 0.8

% I will utilise the corr matrix and distance matrix 
% put them together 
% then index only the ones where corr >= 0.8 
% then make a histogram of that based on distance 

elecSpacingCoef = 200;

distAndCorr = [distanceMatrix(:) * elecSpacingCoef, corrMatrix(:)];

% find distances where correlation > 0.8
distAndHighCorr = distAndCorr(distAndCorr(:, 2) >= 0.8, 1);

figure
histogram(distAndHighCorr, 'EdgeColor', 'white')
aesthetics 
ylabel('Number of connections') 
xlabel('Distance (\mum)')
xlim([200 1700])
set(gca,'TickDir','out'); 
set(gcf, 'Position', [100, 100, 500 * 16/9, 500])
set(gca, 'FontSize', 14)

%% Paper figure:  distance distribution plot multiple histograms 

% same set up, but with medium and low correlation as well 
elecSpacingCoef = 200;
distAndCorr = [distanceMatrix(:) * elecSpacingCoef, corrMatrix(:)];
% find high correlation distacne: > 0.8
distAndHighCorr = distAndCorr(distAndCorr(:, 2) >= 0.8, 1);
% find medium correlation distances: 0.3 - 0.8 
distAndMedCorr = distAndCorr(distAndCorr(:, 2) >= 0.3 & distAndCorr(:, 2) < 0.8, 1);
% find low correlation distances: < 0.3 
distAndLowCorr = distAndCorr(distAndCorr(:, 2) < 0.3, 1);


% one way is to use bar 
% https://uk.mathworks.com/matlabcentral/answers/288261-how-to-get-multiple-groups-plotted-with-histogram
edges = 200:200:1800;

h1 = histcounts(distAndHighCorr, edges) / size(distAndCorr, 1) * 100; % make it a proportion value  
h2 = histcounts(distAndMedCorr, edges) / size(distAndCorr, 1) * 100; 
h3 = histcounts(distAndLowCorr, edges) / size(distAndCorr, 1) * 100; 

figure
bar(edges(1:end-1),[h1; h2; h3]', 'EdgeColor','white')
aesthetics 
ylabel('Proportion of connections (%)') 
xlabel('Distance (\mum)')
% xlim([200 1700])
legend('High correlation (> 0.8)', 'Medium correlation (0.3 - 0.8)', 'Low correlation (< 0.3)')
legend boxoff
set(gca,'TickDir','out'); 
set(gcf, 'Position', [100, 100, 500 * 16/9, 500])
set(gca, 'FontSize', 14)
colormap(viridis)
print(gcf, '-opengl','-depsc', '-r600', '/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/testHist.eps')

% another way is to make transparent histograms
% http://desk.stinkpot.org:8080/tricks/index.php/2006/07/how-to-make-a-transparent-histogram-in-matlab/

figure
histogram(distAndHighCorr, 'EdgeColor', 'white')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75);
hold on
histogram(distAndMedCorr, 'EdgeColor', 'white')
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
aesthetics
legend('High', 'Medium')
legend boxoff

%% Paper figure: distance distribution plot with trendline

% the end product to aim for here is to have just three lines with the
% histfit, and also the data points to show the fit 


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

% same set up, but with medium and low correlation as well 
elecSpacingCoef = 200;
distAndCorr = [distanceMatrix(:) * elecSpacingCoef, corrMatrix(:)];
% find high correlation distacne: > 0.8
distAndHighCorr = distAndCorr(distAndCorr(:, 2) >= 0.8, 1);
% find medium correlation distances: 0.3 - 0.8 
distAndMedCorr = distAndCorr(distAndCorr(:, 2) >= 0.3 & distAndCorr(:, 2) < 0.8, 1);
% find low correlation distances: < 0.3 
distAndLowCorr = distAndCorr(distAndCorr(:, 2) < 0.3, 1);

edges = 200:200:1800;
distAndHighCorrCount = histcounts(distAndHighCorr, edges); 
distAndMedCorrCount = histcounts(distAndMedCorr, edges); 
distAndLowCorrCount = histcounts(distAndLowCorr, edges); 

distAndHighCorrProp = histcounts(distAndHighCorr, edges) / size(distAndCorr, 1) * 100; % make it a proportion value  
distAndMedCorrProp = histcounts(distAndMedCorr, edges) / size(distAndCorr, 1) * 100; 
distAndLowCorrProp = histcounts(distAndLowCorr, edges) / size(distAndCorr, 1) * 100; 


% the actual figure 
figure 
numbins = 8; 
fitmethod = 'gamma'; 

% Low Correlation
h3 = histfit(distAndLowCorr, numbins, fitmethod);
h3Color = [166, 206, 227] / 255;
set(h3(2), 'color', h3Color); 
delete(h3(1))
hold on

% Medium Correlation 
h2 = histfit(distAndMedCorr, numbins, fitmethod); 
h2Color = [178, 223,138] / 255;
set(h2(2),'color', h2Color);
delete(h2(1))
hold on 


% High Correlation 
h1 = histfit(distAndHighCorr, numbins, fitmethod);
h1Color = [31, 120, 180] / 255;
set(h1(2),'color',h1Color);
% remove the bars 
delete(h1(1))
hold on 


% labels and legends 
ylabel('Number of connections') 
xlabel('Distance (\mum)')
aesthetics 
% legend({'High correlation (> 0.8)', 'Medium correlation (0.3 - 0.8)', 'Low correlation (< 0.3)'})
% legend boxoff

% overlay scatter 
dotSize = 50;
scatter(edges(1:8), distAndHighCorrCount, dotSize, h1Color, 'filled')
scatter(edges(1:8), distAndMedCorrCount, dotSize, h2Color, 'filled')
scatter(edges(1:8), distAndLowCorrCount, dotSize, h3Color, 'filled')

xlim([100 1700])
edges = 200:200:1800;
xticks(edges)

% convert from raw count to proportion
% yt = get(gca, 'YTick');
% set(gca, 'YTick', yt, 'YTickLabel', round(yt/size(distAndCorr, 1) * 100 ) )
% ylabel('Proportion of connections (%)')

% only label lines and not the dots 
h = findobj(gca,'Type','line');
legend(h, {'High correlation (> 0.8)', 'Medium correlation (0.3 - 0.8)', 'Low correlation (< 0.3)'})
% legend(h, {'Low correlation (< 0.3)', 'Medium correlation (0.3 - 0.8)', 'High correlation (> 0.8)'})
legend boxoff
 
lineThickness(3); 

set(gca,'TickDir','out'); 
set(gcf, 'Position', [100, 100, 500 * 16/9, 500])
set(gca, 'FontSize', 14)

%% save directly 
print(gcf, '-opengl','-depsc', '-r600', '/media/timothysit/Seagate Expansion Drive1/The_Organoid_Project/figures/paper_figures_first_draft/figure_4/correlation_0503_slice_1_record_2/new_distribution_fits/distribution_gammalfit_high_med_low_slice1_record2_0503_prop.eps')





