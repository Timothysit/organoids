%% Load some data to compare algorithms 

% usually 1209 will do

load('/media/timothysit/Seagate Expansion Drive/The_Mecp2_Project/recordings/mat_files/processed1209/KO_12_09_17-6A_DIV22.mat')

%% raw plot 
subplot(100, 1, [1 25]) 

plot(data)
title('Raw data (1209 DIV22 6A Electrode 58)')
aesthetics()
removeAxis() 

% there is a more elegant way to do this with a for loop 
% and automaticly dividing subplots

%% Prez's method 
[spikeTrain, finalData, threshold] = detectSpikes(data, 'Prez', 4);

subplot(100, 1, [26 30]) % raster 
singleRastPlot(spikeTrain) 
numSpike = sum(spikeTrain);
s = ' ';
title(['Prez: Elliptical, median - 4SD, 1.5ms RP,' s num2str(numSpike) s 'spikes'])
subplot(100, 1, [31 50]) 
plot(finalData); 
% threshold 
hold on; 
plot([1 length(data)], [threshold threshold], '-')
aesthetics()
removeAxis()

%% Manuel's method 
[spikeTrain, finalData, threshold] = detectSpikes(data, 'Manuel', 5);
subplot(100, 1, [51 55]) 
singleRastPlot(spikeTrain) 
numSpike = sum(spikeTrain);
title(['Manuel: Butterworth, mean - 5SD, 2.0ms RP,' s num2str(numSpike) s 'spikes'])
subplot(100, 1, [56 75]) 
plot(finalData); 
% threshold 
hold on; 
plot([1 length(data)], [threshold threshold], '-')
aesthetics()
removeAxis()


%% Tim's method 
[spikeTrain, finalData, threshold] = detectSpikes(data, 'Tim', 8);
subplot(100, 1, [76 80])
singleRastPlot(spikeTrain) 
numSpike = sum(spikeTrain);
title(['Tim: Butterworth, nonlinear energy operator, mean + 12SD, 2.0ms RP,' ... 
    s num2str(numSpike) s 'spikes'])
subplot(100, 1, [81 100]) 
plot(finalData); 
% threshold 
hold on; 
plot([1 length(data)], [threshold threshold], '-')
aesthetics()
removeAxis()


%% Compare Spike V2 



% compare GMM method with threshold-based methods 

% eg. 1209 6A DIV14 e37 

electrode = 37; 

%% Load data 

data = electrodeMatrix(:, electrode); 

%% Number of spikes

spikeIndex = find(pSpike(:, 2) > 0.95); % note p here is for prob

numSpike = [sum(mSpikes(:, electrode)), sum(pSpikes(:, electrode)), ... 
    sum(tSpikes(:, electrode)), length(spikeIndex)];

algorithm = categorical({'Manuel', 'WaveClus', 'NEO', 'GMM'}); 
bar(algorithm, numSpike)

%% Plot the detected waveforms 

fs = 25000;
durationInSec = 2.5 * 10^-3; % the time window to plot trace with spike

figure 
% Manuel spike waveforms 
subplot(1, 5, 1) 
multiplier = 5;
[spikeTrain, finalData, threshold] = detectSpikes(data, 'Manuel', multiplier); 
% [Mspikes, MaverageSpikes] = spikeAlignment(finalData, spikeTrain); 
[Mspikes, MaverageSpikes] = spikeAlignment(finalData, spikeTrain, fs, durationInSec); 
% refractory period is 2.0ms, 50 frames for 25kHz
plotSpikeAlignment(Mspikes);
title(['Manuel RCT paper: ' num2str(sum(spikeTrain)) ' spikes'])
aesthetics

% WaveClus spike waveforms
subplot(1, 5, 2)
multiplier = 4; 
[spikeTrain, finalData, threshold] = detectSpikes(data, 'Prez', multiplier); 
[Pspikes, PaverageSpikes] = spikeAlignment(finalData, spikeTrain); 
% [Pspikes, PaverageSpikes] = spikeAlignment(data, spikeTrain); 
plotSpikeAlignment(Pspikes);
title(['WaveClus: ' num2str(sum(spikeTrain)) ' spikes'])
aesthetics

% NEO Spike waveforms
subplot(1, 5, 3) 
multiplier = 12; 
[spikeTrain, finalData, threshold] = detectSpikes(data, 'Tim', multiplier); 
% [Tspikes, TaverageSpikes] = spikeAlignment(finalData, spikeTrain); 
[Tspikes, TaverageSpikes] = spikeAlignment(data, spikeTrain);
plotSpikeAlignment(Tspikes);
title(['NEO: ' num2str(sum(spikeTrain)) ' spikes'])
aesthetics


% % GMM spike waveform
% subplot(1, 5, 4)
% % spikeIndex = find(spikeTrain == 1); 
% % plot([X(spikeIndex, :)]')
% spikeTrain = gmmDetect(data, 25, 1); 
% [Gspikes, GaverageSpikes] = spikeAlignment(data, spikeTrain); 
% plotSpikeAlignment(Gspikes); 
% title(['GMM: ' num2str(sum(spikeTrain)) ' spikes'])
% aesthetics

% Continuous wavelet transform (default) 
subplot(1, 5, 5) 
% spikeTimes = detect_spikes_wavelet(finalData, 25, [0.5 1.0], 3, 'c', 0.01, 'bior1.5', 0, 0);
% spikeTrain = zeros(size(finalData)); 
% spikeTrain(spikeTimes) = 1;

[spikeTrain, finalData, threshold] = detectSpikes(data, 'cwt', 0);


[cwtSpikes, ] = spikeAlignment(finalData, spikeTrain, fs, durationInSec); 
plotSpikeAlignment(cwtSpikes); 
title(['CWT: ' num2str(sum(spikeTrain)) ' spikes'])
aesthetics 

% I noticed that some of the traces seem to have a rather low amplitude 
% let's focus on those and see if they are spikes 

maxAmp = max(cwtSpikes,[], 2); % get the maximum amplitude of each spike 


figure % look at their distribution
hist(maxAmp)


% seem's like some 'background' like stuff, but can be small spikes!
lowAmpIndex = find(maxAmp < 20); 

figure
plotSpikeAlignment(cwtSpikes(lowAmpIndex, :)); 
aesthtics

% still look quite noisy, must be away to take their max / min to align
% them more properly... 

for putSpike = 1:size(cwtSpikes, 1)
    timeSeries = cwtSpikes(putSpike, :)'; 
    peakIndex = find(max(timeSeries) == timeSeries); 
    if peakIndex > 11 && peakIndex < 40 
        plot(timeSeries(peakIndex-10:peakIndex+10)); 
    end 
end 

%% Compare mSpikes and wSpikes (wavlet) 


[mSpikeTrain, finalData, threshold] = detectSpikes(data, 'Manuel', 5);
[wSpikeTrain, finalData, threshold] = detectSpikes(data, 'cwt', 0);

% look at their total spikes 

bar([sum(MspikeTrain) sum(WspikeTrain)])

% about 1000 more spikes for wavlet algorithm 


% find instances where wavelet found spikes, but m did not
noMatchSpike = find(mSpikeTrain < wSpikeTrain);
% there's actually 2954... most likely because the exact time point found
% is slight different, we can try to downsample it and do this again 

downSize = 50; % let's look at every 50 frames, ie. 2ms
% downMspikeTrain = 


% exclude: 3, 4, 


%% alternatively, compare them over a shorter timewindow 
figure
fs = 25000;
durationInSec = 2.5 * 10^-3; % this is for the spike window

% get the spikes and filtered data 
[mSpikeTrain, mfinalData, threshold] = detectSpikes(data, 'Manuel', 5);
[wSpikeTrain, wfinalData, threshold] = detectSpikes(data, 'cwt', 0, 0);

% eg. look only at the first 5 seconds
signalFrame =  fs * 720; % look at first five seconds 

subplot(2, 2, [1 2])
plot(wfinalData(1:signalFrame), 'HandleVisibility','off'); % don't make legend for this
hold on 
mSpikeLoc = find(mSpikeTrain(1:signalFrame) == 1);
wSpikeLoc = find(wSpikeTrain(1:signalFrame) == 1);

if ~ isempty(mSpikeLoc) % there may be no spikes
    Thrmethod = plot(find(mSpikeTrain(1:signalFrame) == 1), max(max(wfinalData)) + 10, 'Marker','v','MarkerFaceColor','red', 'MarkerSize', 7, 'HandleVisibility','off'); 
end 

if ~ isempty(wSpikeLoc)
    Wavmethod = plot(find(wSpikeTrain(1:signalFrame) == 1), max(max(wfinalData)) + 2, 'Marker','v','MarkerFaceColor','blue', 'MarkerSize', 7, 'HandleVisibility','off');
end
% okay, this is a hard-coded method to get the legend to work
% properly 
% I want only the legend for Thrmethod and Wavmethod, but currently it
% seems to be treating each point as independent, so I can't make a legend
% out of that properly
% 
% hold on 
% plot(mSpikeLoc(1), 14, 'Marker','v','MarkerFaceColor','red', 'MarkerSize', 7, 'DisplayName', 'Threshold method')
% plot(wSpikeLoc(1), 12, 'Marker','v','MarkerFaceColor','blue', 'MarkerSize', 7, 'DisplayName', 'Wavelet method')
% 
% legend show 
% % lg = findobj(gcf, 'Type', 'Legend'); 
% % lineInLegend = findobj(lg,'type','line','linestyle','-');
% % set(lineInLegend,'visible','off');
% legend boxoff

aesthetics

title('1209 6A DIV 22 e37, 0 - 720s. Butterworth filtered')

% will be nice to have arrows on top of this plot to indicate the 
% time when each algorithm detected spikes

subplot(2, 2, 3) 
[mSpikes, ] = spikeAlignment(mfinalData(1:signalFrame), mSpikeTrain(1:signalFrame), fs, durationInSec); 
plotSpikeAlignment(mSpikes, 'peak'); 
aesthetics
title(['Spike traces: Threshold method. ' num2str(length(mSpikeLoc)), ' spikes'])

subplot(2, 2, 4) 
[cwtSpikes, ] = spikeAlignment(wfinalData(1:signalFrame), wSpikeTrain(1:signalFrame), fs, durationInSec); 
plotSpikeAlignment(cwtSpikes, 'peak'); 
aesthetics
title(['Spike traces: Wavelet method. ' num2str(length(wSpikeLoc)), ' spikes'])

%% 1209 4A DIV 14 case study


load('/media/timothysit/Seagate Expansion Drive/The_Mecp2_Project/recordings/mat_files/goodFiles/HT_12_09_17-4A_DIV14.mat')



% 23 - defeinitely active 
% 22 - putative low activity 
% 4 - definitely inactive --> grounded 

%% Let's tune the L parameter in the CWT function 

% sorry this is quite messy and very hard-coded
figure
fs = 25000;

% eg. look only at the first 5 seconds
signalFramRange =  (360 * fs):(361 *fs); % 5 seconds middle of the recording 

subplot(3, 1, 1) 
L = -0.2; 
[wSpikeTrain, wfinalData, threshold] = detectSpikes(data, 'cwt', 0, L);
plot(wfinalData(signalFramRange), 'HandleVisibility','off'); % don't make legend for this
hold on 
wSpikeLoc = find(wSpikeTrain(signalFramRange) == 1);

if ~ isempty(wSpikeLoc)
    Wavmethod = plot(find(wSpikeTrain(signalFramRange) == 1), max(max(wfinalData)) + 2, 'Marker','v','MarkerFaceColor','blue', 'MarkerSize', 7, 'HandleVisibility','off');
end
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3, 1, 2) 

L = 0; 
[wSpikeTrain, wfinalData, threshold] = detectSpikes(data, 'cwt', 0, L);
plot(wfinalData(signalFramRange), 'HandleVisibility','off'); % don't make legend for this
hold on 
wSpikeLoc = find(wSpikeTrain(signalFramRange) == 1);

if ~ isempty(wSpikeLoc)
    Wavmethod = plot(find(wSpikeTrain(signalFramRange) == 1), max(max(wfinalData)) + 2, 'Marker','v','MarkerFaceColor','blue', 'MarkerSize', 7, 'HandleVisibility','off');
end

%%%%%%%%%%%%%%%%%%%%%
subplot(3, 1, 3)

L = 0.2; 
[wSpikeTrain, wfinalData, threshold] = detectSpikes(data, 'cwt', 0, L);
plot(wfinalData(signalFramRange), 'HandleVisibility','off'); % don't make legend for this
hold on 
wSpikeLoc = find(wSpikeTrain(signalFramRange) == 1);

if ~ isempty(wSpikeLoc)
    Wavmethod = plot(find(wSpikeTrain(signalFramRange) == 1), max(max(wfinalData)) + 2, 'Marker','v','MarkerFaceColor','blue', 'MarkerSize', 7, 'HandleVisibility','off');
end

%% Plot parameter against number of spikes detected 

% for threshold method; plot threshold against number of spikes detected 
% for wavelet method; plot L (loss ratio) against number of spikes detected

% plan is to do this for differnet types of channels 

% 1209 4A DIV 14
% grounded electrode: 4
% active electrode, with cell: 23 
% pseudo-active electroed, without cell: 34 (or 35, 36)
% possibly active, with cell: 17, 31 

% 1209 6A DIV 22 
% grounded electrode 4
% very very active, with large cluster: 37



data = electrodeMatrix(:, 37); 

subplot(1, 2, 1)


% detectionMethods = {'Prez', 'Manuel', 'Tim'}; 
detectionMethods = {'Prez', 'Manuel'}; 
multiplier = linspace(4, 7, 9); 
spikeStore = zeros(1, length(multiplier)); 

for method = 1:length(detectionMethods)
    for multIndex = 1:length(multiplier)
        % feed the threshold into detection algorithm
        [spikeTrain, finalData, threshold] = ... 
        detectSpikes(data, detectionMethods{method}, multiplier(multIndex)); % for one electrode 
        % test for all electrodes (TODO)
        numSpikes = sum(spikeTrain);
        spikeStore(multIndex) = numSpikes;
    end
    plot(multiplier, spikeStore)
    % plot(log10(spikeStore))
    % plot(gradient(spikeStore)); 
    % there is suggestion that finding the peak of the gradient 
    % is a good way to optimise the threshold value
    hold on
end
% title('Threshold-based methods')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlabel('Threshold multiplier')
% ylabel('Log(Number of spikes)')
ylabel('Number of detected spikes')
aesthetics()
% legend('Elliptical', 'Butterworth', 'Butterworth + NEO')
legend('Elliptical-median', 'Butterworth-mean')
legend boxoff  

% set(gcf, 'position', [100 100 500*1.6 500])

subplot(1, 2, 2) 

Lparam = linspace(-0.2, 0.2, 9);
spikeStore = zeros(1, length(Lparam)); 
for Lindex = 1:length(Lparam)
    L = Lparam(Lindex);
    [spikeTrain, finalData, threshold] = detectSpikes(data, 'cwt', 0, L);
    spikeStore(Lindex) = sum(spikeTrain);
end 

plot(Lparam, spikeStore)
aesthetics()
lineThickness(2)
xlabel('Loss ratio') 
ylabel('Number of detected spikes') 
legend('Wavelet')
legend boxoff
% title('Wavelet based method')
% set(gcf, 'position', [100 100 500*1.6 500])

suptitle('1209 4A DIV14 e4 (grounded)')
set(gcf, 'Position', [100, 100, 650, 500])


%% Just plot the filtered signal and the detected spike times (as triangles) 
% Select data to load here.
% uiopen('/media/timothysit/Seagate Expansion Drive/The_Mecp2_Project/recordings/mat_files/')


figure
fs = 25000;
durationInSec = 2.5 * 10^-3; % this is for the spike window

% get the spikes and filtered data 
[mSpikeTrain, mfinalData, threshold] = detectSpikes(data, 'Manuel', 5, 0);
[wSpikeTrain, wfinalData, threshold] = detectSpikes(data, 'cwt', 0, -0.1254);

% specify time range to look at
% signalFrame =  fs * 720; 
signalFrameRange = (60 * fs):(90 * fs); % look at entire recording
% signalFrameRange = [360 * fs: 365 * fs];  

plot(wfinalData(signalFrameRange), 'HandleVisibility','off'); % don't make legend for this
hold on 
mSpikeLoc = find(mSpikeTrain(signalFrameRange) == 1);
wSpikeLoc = find(wSpikeTrain(signalFrameRange) == 1);

if ~ isempty(mSpikeLoc) % there may be no spikes
    Thrmethod = plot(find(mSpikeTrain(signalFrameRange) == 1), max(max(wfinalData)) + 10, 'Marker','v','MarkerFaceColor','red', 'MarkerEdgeColor', 'none', 'MarkerSize',  7, 'HandleVisibility','off'); 
end 

if ~ isempty(wSpikeLoc)
    Wavmethod = plot(find(wSpikeTrain(signalFrameRange) == 1), max(max(wfinalData)) + 2, 'Marker','v','MarkerFaceColor','blue', 'MarkerEdgeColor', 'none', 'MarkerSize', 7, 'HandleVisibility','off');
end
% okay, this is a hard-coded method to get the legend to work
% properly 
% I want only the legend for Thrmethod and Wavmethod, but currently it
% seems to be treating each point as independent, so I can't make a legend
% out of that properly

hold on 
plot(mSpikeLoc(1), max(max(wfinalData)) + 10, 'Marker','v','MarkerFaceColor','red', 'MarkerEdgeColor', 'none', 'MarkerSize', 7, 'DisplayName', 'Threshold method')
plot(wSpikeLoc(1), max(max(wfinalData)) + 2, 'Marker','v','MarkerFaceColor','blue', 'MarkerEdgeColor', 'none', 'MarkerSize', 7, 'DisplayName', 'Wavelet method')

legend('location', 'northeastoutside')
legend show 
lg = findobj(gcf, 'Type', 'Legend'); 
lineInLegend = findobj(lg,'type','line','linestyle','-');
set(lineInLegend,'visible','off');
legend boxoff

% xlim([1 length(signalFrameRange)]) % tight x-axis

aesthetics
xlabel('Frames (25kHz)') 
ylabel('Amplitude')
% title('0411 5B DIV 28 e58, 0 - 720s. Butterworth filtered')
height = 500; 
% width = height * 1.618; % golden ratio
width = height * 3; % for long display
set(gcf, 'Position', [100 100 width height]);

% some embellishment 
removeAxis 
% change legend size
lg.FontSize = 14;
lg.Location = 'southwest';
sb = scalebar;
sb.Position = [0 -10];
sb.YUnit = 'arb. units';
sb.XUnit = '2 seconds';
sb.YLen = 4;

%% Plot the detected spike traces: Threshold method 
figure
[mSpikes, ] = spikeAlignment(mfinalData(signalFrameRange), mSpikeTrain(signalFrameRange), fs, durationInSec); 
% plotSpikeAlignment(mSpikes, 'peak'); 
plotSpikeAlignment(mSpikes, 'simple', 25000, 1.5 * 10^-3); 
xlabel('Frames (25kHz)')
ylabel('Amplitude')
aesthetics
set(gcf, 'Position', [100 100 width/2 height]);
% title(['Spike traces: Threshold method. ' num2str(length(mSpikeLoc)), ' spikes'])

removeAxis
sb = scalebar;
sb.YUnit = 'arb. units';
sb.XUnit = 'ms';
sb.YLen = 3;
sb.XLen = 10;

%% Plot the detected spike traces: Wavelet method 
figure 
[cwtSpikes, ] = spikeAlignment(wfinalData(signalFrameRange), wSpikeTrain(signalFrameRange), fs, durationInSec); 
plotSpikeAlignment(cwtSpikes, 'peak', 25000, 1.5 * 10^-3); 

xlabel('Frames (25kHz)')
ylabel('Amplitude')
aesthetics
% title(['Spike traces: Wavelet method. ' num2str(length(wSpikeLoc)), ' spikes'])
set(gcf, 'Position', [100 100 width/2 height]);

removeAxis
sb = scalebar;
sb.Position = [0 -8];
sb.YUnit = 'arb. units';
sb.XUnit = 'ms';
sb.YLen = 3;
sb.XLen = 10;
lineThickness(3)

%% Plot the number of detected spikes from each method - one electrode 


%% Get cell of m5Spikes
dirName = '/media/timothysit/Seagate Expansion Drive/The_Mecp2_Project/recordings/mat_files/goodFiles_DIV14_DIV35/spikeFiles/m5Spikes/';
files = dir(dirName);
files = files(~ismember({files.name},{'.','..'}));

for file = 1:length(files) 
   load([dirName files(file).name]) 
   m5SpikeElectrode = sum(m5Spikes);
   mSpikes{file} = m5SpikeElectrode;
end

%% Get cell of tSpikes of pSpikes
% TODO: combine all of them into one big nested for loop will give a more
% elegand solution.

dirName = '/media/timothysit/Seagate Expansion Drive/The_Mecp2_Project/recordings/mat_files/goodFiles_DIV14_DIV35/spikeFiles/p_and_t_spikes/';
files = dir(dirName);
files = files(~ismember({files.name},{'.','..'}));

for file = 1:length(files) 
   load([dirName files(file).name]) 
   t12SpikeElectrode = sum(t12Spikes);
   tSpikes{file} = t12SpikeElectrode; 
   
   p4SpikeElectrode = sum(p4Spikes);
   pSpikes{file} = p4SpikeElectrode; 
end

%% wavelet number spike values 

dirName = '/media/timothysit/Seagate Expansion Drive/The_Mecp2_Project/recordings/mat_files/goodFiles_DIV14_DIV35/spikeFiles/cwt_neg01254_noRefrac/'; 
files = dir(dirName);
files = files(~ismember({files.name},{'.','..','featMatrix'}));
for file = 1:length(files) 
   load([dirName files(file).name]) 
   cwtSpikeElectrode = sum(cNeg01254_Spikes);
   cwtSpikes{file} = cwtSpikeElectrode; 
end

%% Combine all of them into vectors so that I can find the mean and SEM

% requires barwitherror from barChartUtil folder

mSpikes = cell2mat(mSpikes); 
tSpikes = cell2mat(tSpikes); 
pSpikes = cell2mat(pSpikes); 
cwtSpikes = cell2mat(cwtSpikes); 



spikeCountmeans = [mean(pSpikes), mean(mSpikes), mean(tSpikes), mean(cwtSpikes)]; % mean
getSem = @(data) std(data)/sqrt(length(data)); % standard error of mean
spikeCountSEM = [getSem(pSpikes), getSem(mSpikes), getSem(tSpikes), getSem(cwtSpikes)]; 




%%  plot bar graph
figure 
% bar(spikeCountmeans); 
methods = {'A','B','C', 'D'}; 
b = barwitherr(spikeCountSEM, spikeCountmeans, 'FaceColor',[0 .5 .5], 'EdgeColor',[1 1 1]);
aesthetics 
set(gca,'XTickLabel',methods)
set(gca, 'TickDir', 'out') % tickmarks point outwards
ylabel('Mean numebr of spikes') 
xlabel('Spike detection method')
set(gca,'FontSize',14)
set(gcf, 'Position', [100 100 900 800]);

% remodify colour of bars to match that in spikeDataParamTune 
% I have decided to use the tetradic color harmony (rectangle) scheme
% See: https://luminous-landscape.com/color-harmonies-4-cool-warm-split-tetradic-and-square/
colours = {'red', 'purple', 'yellow green', 'blue green'};



% A: Prez method
% B: Manuel method 
% C: NEO method 
% D: Wavlet method c = -0.1254 

% should also include the wavelet method 

%% Plot each bar individually so I can give them different colours 
methods = {'A','B','C', 'D'}; 
figure 
hold on
colours = {'red', 'purple', 'yellow green', 'blue green'};
for i = 1:length(methods)
    % b = barwitherr(spikeCountSEM(i), spikeCountmeans(i), 'FaceColor',[0 0 0], 'EdgeColor',[1 1 1]);
    % hold on
    h = bar(i, spikeCountmeans(i));
    set(h, 'FaceColor', rgb(colours{i}), 'EdgeColor', [1 1 1])
    % h=bar(i,mydata(i));
    hold on
end
aesthetics

% error bar 
errorbar(1:length(spikeCountmeans),spikeCountmeans,spikeCountSEM,'.', 'LineWidth', 2, 'Color', 'black')
xticks([1 2 3 4])
set(gca,'XTickLabel',methods)
set(gca, 'TickDir', 'out') % tickmarks point outwards
ylabel('Mean numebr of spikes') 
xlabel('Spike detection method')
set(gca,'FontSize',14)
% set(gcf, 'Position', [100 100 900 800]);
set(gcf, 'Position', [100, 100, 650, 500])
%% performance on electrode 4 only (grounded electrode)

% this may not work... since some of the electrodes have electrode 4
% removed, oh wait.. yes, I think I removed all of them from the current
% spike files...
mGrounded = 0; 
pGrounded = 0; 
tGrounded = 0; 
cwtGrounded = 0; 
for mea = 1:length(mSpikes)
    mGrounded = mGrounded + sum(mSpikes{mea}(:, 4)); 
    pGrounded = pGrounded + sum(pSpikes{mea}(:, 4)); 
    tGrounded = tGrounded + sum(tSpikes{mea}(:, 4)); 
    cwtGrounded = cwtGrounded + sum(cwtSpikes{mea}(:, 4)); 
end 

% plot 
figure 
methods = categorical({'A','B','C'});
% barwitherr(spikeCountSEM, spikeCountmeans, 'FaceColor',[0 .5 .5], 'EdgeColor',[1 1 1]);
aesthetics 


