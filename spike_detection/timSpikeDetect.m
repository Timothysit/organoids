% Tim's spike detection method 

figure 

%% Raw data 
subplot(100, 1, [1 40])
plot(data)
title('Raw data')
aesthetics()

%% butterworth filter
lowpass = 600; 
highpass = 8000; 
fs = 25000; 
wn = [lowpass highpass] / (fs / 2); 
filterOrder = 3;
[b, a] = butter(filterOrder, wn); 
filteredData = filtfilt(b, a, double(data)); 

%% NEO, threshold and spikes

% NEO by calling snle
y_snle = snle(filteredData', 1); 
m = mean(y_snle); 
s = std(y_snle); 

% 20171213: in accordance to the original paper (Mukhopadhyay and Ray
% 1998), I will use a scaled mean rather than a standard dev based
% threshold.

% some suggestion that 8 is the value to use here (Gibson et al 2008, but
% will be dependent of the data you have). Use compareSpikeDetect.m and
% plotSpikeAlignment.m to verify that the detected spikes do look like
% spikes

multiplier = 8; % this is the crux of the detection 
% threshold = m + multiplier*s; 
threshold = m * multiplier;
spikeTrain = y_snle > threshold; 
% this is a much large std than what others had to use...
% but this is because we used NEO
spikeTrain = double(spikeTrain);
numSpike = sum(spikeTrain);


% plot NEO, then raster plot 

% raster plot 
subplot(100, 1, [48 51])
singleRastPlot(spikeTrain') 

subplot(100, 1, [56 100])
plot(y_snle)
% titleText = '';
title('Butteworth, Non-linear Energy Operator, 12 SD')
hold on 
plot([1 length(data)], [threshold threshold], '-')
aesthetics()







