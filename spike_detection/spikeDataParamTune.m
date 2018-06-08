% detection algorithm parameter tuning 

% Author: Tim Sit 

% one way to tune the spike threshold is to look at the number of spikes
% detected in a grounded electrode. And set the threshold so that no spikes
% are detected from the grounded electrode (this will of course be a
% conservative method)

%% Load data 

% load grounded electrode 
% let's try 1209 6A DIV22 

load('/media/timothysit/Seagate Expansion Drive/The_Mecp2_Project/recordings/mat_files/goodFiles/KO_12_09_17-6A_DIV22.mat')

data = electrodeMatrix(:, 4); % look at grounded electrode.

%% plot the waveform of the grounded electrode

figure 
gridTrace(electrodeMatrix, 1000);
ylim([-7804, -7792])

figure 
plot(data) 
aesthetics
% removeAxis


%% Threshold and number of spikes 

spikeStore = zeros(1, 10); 
detectionMethods = {'Prez', 'Manuel', 'Tim'}; 
colours = {'red', 'purple', 'yellow green'};
for method = 1:length(detectionMethods)
    for multiplier = 1:10
    % feed the threshold into detection algorithm
        [spikeTrain, finalData, threshold] = ... 
            detectSpikes(data, detectionMethods{method}, multiplier); % for one electrode 
        % test for all electrodes (TODO)
        numSpikes = sum(spikeTrain);
        spikeStore(multiplier) = numSpikes;
    end
    % plot(spikeStore)
    plot(spikeStore, 'color', rgb(colours{method})) % require xkcdRGB
    % plot(log10(spikeStore))
    % plot(gradient(spikeStore)); 
    % there is suggestion that finding the peak of the gradient 
    % is a good way to optimise the threshold value
    hold on
end
% title('2309 4F DIV 27 E1')
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlabel('Threshold multiplier')
% ylabel('Log(Number of spikes)')
ylabel('Number of detected spikes')
lineThickness(2.5)
aesthetics()
legend('A', 'B', 'C')
legend boxoff  
set(gca,'TickDir','out');
set(gca,'FontSize', 13)
set(gcf, 'Position', [100, 100, 650, 500])

%% Wavelet method parameter

costParam = [-0.188, -0.1254, 0, 0.1254, 0.188];
spikeStore = zeros(1, length(costParam));
multiplier = 0; % wavelet method doens't depend on the multiplier param
figure 
 for cc = 1:length(costParam)
    % feed the threshold into detection algorithm
        [spikeTrain, finalData, threshold] = ... 
            detectSpikes(data, 'cwt', multiplier, cc); % for one electrode 
        % test for all electrodes (TODO)
        numSpikes = sum(spikeTrain);
        spikeStore(cc) = numSpikes;
 end
% scatter(costParam, spikeStore) 
% hold on
plot(costParam, spikeStore, 'color', rgb('blue green'))
aesthetics
legend('D') 
legend boxoff
xticks(costParam)
xlabel('Cost parameter')
ylabel('Number of detected spikes')
lineThickness(2.5)
set(gca,'TickDir','out');
set(gca,'FontSize', 13)
set(gcf, 'Position', [100, 100, 650, 500])


% plot(log10(spikeStore))
% plot(gradient(spikeStore)); 
% there is suggestion that finding the peak of the gradient 
% is a good way to optimise the threshold value




%% Mean vs Median 


%% Refractory period 
% for refPeriod = 0:50 
%             [spikeTrain, finalData, threshold] = ... 
%             detectSpikes(data, detectionMethods{method}, multiplier);
%         numSpikes = sum(spikeTrain);
%         spikeStore(multiplier) = numSpikes;
% end 

%% 