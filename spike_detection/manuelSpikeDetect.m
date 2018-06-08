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

% finding threshold and spikes
m = mean(filteredData); 
s = std(filteredData); 
multiplier = 5;
threshold = m - multiplier*s; 
spikeTrain = filteredData < threshold; % negative threshold

% impose refractory period
refPeriod = 2.0 * 10^-3 * fs; % 2ms 
% I think there is a more efficient/elegant way to do this, but I haven't 
% taken time to think about it yet 
spikeTrain = double(spikeTrain);
for i = 1:length(spikeTrain)
    if spikeTrain(i) == 1 
        refStart = i + 1; % start of refractory period 
        refEnd = round(i + refPeriod); % end of refractory period
        spikeTrain(refStart:refEnd) = 0; 
    end 
end 


numSpike = sum(spikeTrain);

% raster plot 
subplot(100, 1, [48 51])
singleRastPlot(spikeTrain) 

subplot(100, 1, [56 100])
plot(filteredData); 
title('4th order Elliptical Filter 300 - 8000Hz, 5SD threshold')
% threshold 
hold on; 
plot([1 length(data)], [threshold threshold], '-')
aesthetics()