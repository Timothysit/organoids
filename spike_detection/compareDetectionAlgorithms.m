% compare detection techniques 


figure

%% Raw plot
% subplot(1, 100, [1 25])


%% Prez's method 
% subplot(1, 100, [26 30]) % raster 
% subplot(1, 100, [31 50]) % threshold


%% Manuell's method 
% subplot(1, 100, [51, 55]) % raster 
% subplot(1, 100, [56 - 75]) % threshold 


%% Tim's method 
% subplot(1, 100, [76 - 80]) % raster 
% subplot(1, 100, [81 - 100]) threshold


%% Compare Detected spikes waveforms
% load the electrode matrix file you want 
% data = electrodeMAtrix(:, electrode #)


figure 
% Manuel spike waveforms 
subplot(3, 2, 1) 
multiplier = 5;
[spikeTrain, finalData, threshold] = detectSpikes(data, 'Manuel', multiplier); 
[Mspikes, MaverageSpikes] = spikeAlignment(data, spikeTrain); 
plotSpikeAlignment(Mspikes);
hold off 

% Prez spike waveforms
subplot(3, 2, 3)
multiplier = 4; 
[spikeTrain, finalData, threshold] = detectSpikes(data, 'Prez', multiplier); 
[Pspikes, PaverageSpikes] = spikeAlignment(data, spikeTrain); 
plotSpikeAlignment(Pspikes);


% Tim Spike waveforms 
subplot(3, 2, 5) 
multiplier = 12; 
[spikeTrain, finalData, threshold] = detectSpikes(data, 'Tim', multiplier); 
[Tspikes, TaverageSpikes] = spikeAlignment(data, spikeTrain); 
plotSpikeAlignment(Tspikes);

% everyone's average spike waveform 
subplot(3, 2, [2, 4, 6])
title('Average spike form (1209 6A DIV 22 E 11)')
window = -floor(length(MaverageSpikes)/2):floor(length(MaverageSpikes)/2);
plot(window, MaverageSpikes)
hold on 
plot(window, PaverageSpikes) 
hold on 
plot(window, TaverageSpikes) 
legend('Manuel (1585 spikes, 5SD)', 'Prez (10932 spikes, 4SD)', 'Tim (4223 spikes, 12SD)')
legend boxoff 
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
aesthetics()

%% Average waveform with varying threshold 

figure 
for multiplier = 3:10
    [spikeTrain, finalData, threshold] = detectSpikes(data, 'Tim', multiplier); 
    [spikes, averageSpikes] = spikeAlignment(data, spikeTrain); 
    window = -floor(length(averageSpikes)/2):floor(length(averageSpikes)/2);
    plot(window, averageSpikes) 
    hold on 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    aesthetics()
end 
title('Tim, 2.0ms')
xlabel('Time since spike on set (in frames)')
legend(string(3:10), 'Location', 'southwest')
legend boxoff 

%% Plot the three together (average spike with different thresholds)

subplot(1, 3, 1) 
for multiplier = 3:10
    [spikeTrain, finalData, threshold] = detectSpikes(data, 'Prez', multiplier); 
    [spikes, averageSpikes] = spikeAlignment(data, spikeTrain); 
    window = -floor(length(averageSpikes)/2):floor(length(averageSpikes)/2);
    plot(window, averageSpikes) 
    hold on 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    aesthetics()
    removeYAxis
end 
title('Prez, 1.5ms')
xlabel('Time since spike onset (in frames)')
leg = legend(string(3:10), 'Location', 'southwest'); 
title(leg, 'Threshold multiplier')
legend boxoff


subplot(1, 3, 2) 
for multiplier = 3:10
    [spikeTrain, finalData, threshold] = detectSpikes(data, 'Manuel', multiplier); 
    [spikes, averageSpikes] = spikeAlignment(data, spikeTrain); 
    window = -floor(length(averageSpikes)/2):floor(length(averageSpikes)/2);
    plot(window, averageSpikes) 
    hold on 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    aesthetics()
    removeYAxis
end 
title('Manuel, 2.0ms')
xlabel('Time since spike onset (in frames)')
leg = legend(string(3:10), 'Location', 'southwest'); 
title(leg, 'Threshold multiplier')
legend boxoff



subplot(1, 3, 3) 
for multiplier = 3:10
    [spikeTrain, finalData, threshold] = detectSpikes(data, 'Tim', multiplier); 
    [spikes, averageSpikes] = spikeAlignment(data, spikeTrain); 
    window = -floor(length(averageSpikes)/2):floor(length(averageSpikes)/2);
    plot(window, averageSpikes) 
    hold on 
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    aesthetics()
    removeYAxis
end 
title('Tim, 2.0ms')
xlabel('Time since spike onset (in frames)')
leg = legend(string(3:10), 'Location', 'southwest'); 
title(leg, 'Threshold multiplier')
legend boxoff 







