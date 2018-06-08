function [spikeTrain, finalData, threshold] = detectSpikes(data, method, multiplier, L)

% input 
    % data: a n x 1 vector containing the signal, where n is the number of
    % samples 
    % method: string value specifying the spike detection method to use 
    % multiplier: the threshold multiplier to use for your chosen method 
    % L: loss parameter for wavelet method, won't matter if other methods
    % used. Default is zero as recommended by the creator.

% Author: Tim Sit, sitpakhang@gmail.com 
% Last update: 20180503

    
if ~exist('L')
    L = 0; 
end 
    

    

%% General paramters

fs = 25000; % sampling rate

% artifact removal 

% filteredData = artifactThresh(filteredData,validMask,thresh); 

%% Prez's method 

if strcmp(method,'Prez') 
    par.detect_order = 4; % default, no idea why specifically this number
    par.ref_ms = 1.5; % refractor period in ms, not sure when this is going to be used...
    refPeriod = par.ref_ms * 10^(-3) * fs; % covert to frames
    fmin_detect = 300; 
    fmax_detect = 8000;
    [b,a] = ellip(par.detect_order,0.1,40,[fmin_detect fmax_detect]*2/fs);
    % FiltFiltM does the same thing, but runs slightly faster
    % the parameters are default found in wave_clus
    filteredData = filtfilt(b, a, double(data)); 

    % finding threshold and spikes
    med = median(filteredData); 
    s = std(filteredData); 
    % multiplier = 4; % default value
    threshold = med - multiplier*s; 
    spikeTrain = filteredData < threshold; % negative threshold
    spikeTrain = double(spikeTrain);
    finalData = filteredData;
    
   %  although it wasn't mentioned in the thesis, I think he implemented
   % the default Wave_clus 1.5 ms refractory period
    for i = 1:length(spikeTrain)
       if spikeTrain(i) == 1 
           refStart = i + 1; % start of refractory period 
           refEnd = round(i + refPeriod); % end of refractory period
           if refEnd > length(spikeTrain) 
               % prevents extending the vector
               % in the case there is a spike at the end of the recording
               spikeTrain(refStart:length(spikeTrain)) = 0; 
           else 
                spikeTrain(refStart:refEnd) = 0; 
           end
       end 
    end 
    
    
    % alternative version to speed the above up to avoid unncessary loops
    % The spike count is different... need to spend some time to dissect it
    % Run time is actually quite similar
%     L = length(spikeTrain); 
%     spikeFrames = find(spikeTrain == 1); 
%     for i = 1:length(spikeFrames)
%         if spikeTrain(spikeFrames(i)) == 1 
%             % check that the spike havne't already been removed previously 
%             refStart = spikeFrames(i) + 1; 
%             refEnd = round(spikeFrames(i)  + refPeriod); 
%             spikeTrain(refStart:refEnd) = 0; 
%         end 
%     end 
%     spikeTrain = spikeTrain(1:L); 
    % remove refractory 0s added to end of recording
    % in case there is spike at the end of recording 
% based on this: 
% https://www.mathworks.com/matlabcentral/fileexchange/55227-automatic-objective-neuronal-spike-detection?focused=8345812&tab=function
    
end 


%% implement different methods to detect spikes  
if strcmp(method,'Tim')
    % code inspired by Gaidica 
    % http://gaidi.ca/weblog/extracting-spikes-from-neural-electrophysiology-in-matlab
    % butterworth filter 
    lowpass = 600; 
    highpass = 8000; 
    wn = [lowpass highpass] / (fs / 2); 
    filterOrder = 3;
    [b, a] = butter(filterOrder, wn); 
    filteredData = filtfilt(b, a, double(data)); 
    % NEO by calling snle
    y_snle = snle(filteredData', 1); 
    m = mean(y_snle); 
    s = std(y_snle); 
    % multiplier = 12; % this is the crux of the detection 
    
    % 20171123: I have changed this to match with the original
    % implementation (Mukhodpadhyay and Ray 1998);
    % to use a scaled mean as the threshold rather than a
    % standard deviation based approach
    
    % 20180413: referred back to this just for historical comparison for
    % thesis
    threshold = m + multiplier*s; 
    % threshold = m * multiplier;
    spikeTrain = y_snle > threshold; 
    % this is a much large std than what others had to use...
    % but this is because we used NEO
    spikeTrain = double(spikeTrain)';
    
    % refractory period 
    refPeriod = 2.0 * 10^-3 * fs; % 2ms 
    for i = 1:length(spikeTrain)
       if spikeTrain(i) == 1 
           refStart = i + 1; % start of refractory period 
           refEnd = round(i + refPeriod); % end of refractory period
           if refEnd > length(spikeTrain)
               spikeTrain(refStart:length(spikeTrain)) = 0; 
           else 
               spikeTrain(refStart:refEnd) = 0; 
           end 
       end 
    end     
    
    
    finalData = y_snle;
end 

% Continuous Wavelet Transform 
% http://cbmspc.eng.uci.edu/SOFTWARE/SPIKEDETECTION/tutorial/tutorial.html

%% M.Schroter 2015 Spike Detection Procedure 

% bandpass filter : 3rd order Butterworth, 600 - 8000Hz 
% threshold of 5 x SD below backgroudn noise for each channel 
% impose 2 ms refractory period (Wagennaar et al 2006) 

% burst detection via ...
% spike times downsample to 1KHz 
% activity of all  electrodes averaged over windows of 10ms into one vector
% vector serached for clusters of activity ( <60ms inter-event internval) 
% if activity within cluster occured on at least 6 electrodes and contain
% at least 50 spike, then population spike 
% the numbers seem quite arbitrary to me...
if strcmp(method,'Manuel')
    % butterworth filter 
    lowpass = 600; 
    highpass = 8000; 
    wn = [lowpass highpass] / (fs / 2); 
    filterOrder = 3;
    [b, a] = butter(filterOrder, wn); 
    filteredData = filtfilt(b, a, double(data)); 

    % finding threshold and spikes
    m = mean(filteredData); 
    s = std(filteredData); 
    % multiplier = 5;
    threshold = m - multiplier*s; 
    % negThreshold = m - 8 * s; % maximum threshold, a simple artefact removal method 
    spikeTrain = filteredData < threshold; 
    
   
 

    % impose refractory period
    refPeriod = 2.0 * 10^-3 * fs; % 2ms 
    % I think there is a more efficient/elegant way to do this, but I haven't 
    % taken time to think about it yet 
    spikeTrain = double(spikeTrain);
    finalData = filteredData;
    % refractory period
    for i = 1:length(spikeTrain)
       if spikeTrain(i) == 1 
           refStart = i + 1; % start of refractory period 
           refEnd = round(i + refPeriod); % end of refractory period
           if refEnd > length(spikeTrain)
               spikeTrain(refStart:length(spikeTrain)) = 0; 
           else 
               spikeTrain(refStart:refEnd) = 0; 
           end 
       end 
    end 
end 

%%  Z. Nenadic and J.W. Burdick 2005 
% , Spike detection using the 
%   continuous wavelet transform, IEEE T. Bio-med. Eng., vol. 52,
%   pp. 74-87, 2005.


if strcmp(method,'cwt')
    
    
    % Filter signal 
    lowpass = 600; 
    highpass = 8000; 
    wn = [lowpass highpass] / (fs / 2); 
    filterOrder = 3;
    [b, a] = butter(filterOrder, wn); 
    filteredData = filtfilt(b, a, double(data)); 
    
    % input sampling frequency as kHz
    Wid = [0.5 1.0]; % 1 x 2 vector of expected minimum and maximum width [msec] of transient 
    %  to be detected Wid=[Wmin Wmax]. For most practical purposes Wid=[0.5 1.0];
    
    % TS 20180228, let's try to be slightly more generous... 
    % Wid = [0.2 1.0]; result: same, no change in spike count
    
    Ns = 5; % Ns - (scalar): the number of scales to use in detection (Ns >= 2);
    option = 'c'; % the action taken when no coefficients survive hard thresholding 
    %   'c' means conservative and returns no spikes if P(S) is found to be 0
    %   'l' means assume P(S) as a vague prior (see the original reference)
    % L = -0.2; % L is the factor that multiplies [cost of comission]/[cost of omission].
    %   For most practical purposes -0.2 <= L <= 0.2. Larger L --> omissions
    %   likely, smaller L --> false positives likely. For unsupervised
    %   detection, the suggested value of L is close to 0.  
    
    wname = 'bior1.5'; 
    
    %   wname - (string): the name of wavelet family in use
    %           'bior1.5' - biorthogonal
    %           'bior1.3' - biorthogonal
    %           'db2'     - Daubechies
    %           'sym2'    - symmlet
    %           'haar'    - Haar function
    %   Note: sym2 and db2 differ only by sign --> they produce the same
    %   result;
    
    PltFlg = 0; 
    CmtFlg = 0; 
    %   PltFlg - (integer) is the plot flag: 
    %   PltFlg = 1 --> generate figures, otherwise do not;
    %  
    %   CmtFlg - (integer) is the comment flag, 
    %   CmtFlg = 1 --> display comments, otherwise do not;

    spikeFrames = detect_spikes_wavelet(filteredData, fs/1000, Wid, Ns, option, L, wname, PltFlg, CmtFlg); 
    
    spikeTrain = zeros(size(data)); 
    spikeTrain(spikeFrames) = 1;
    
     % impose refractory period
    % refPeriod = 2.0 * 10^-3 * fs; % 2ms 
    % I think there is a more efficient/elegant way to do this, but I haven't 
    % taken time to think about it yet 
    % refractory period
    % for i = 1:length(spikeTrain)
    %    if spikeTrain(i) == 1 
    %        refStart = i + 1; % start of refractory period 
    %        refEnd = round(i + refPeriod); % end of refractory period
    %        if refEnd > length(spikeTrain)
    %            spikeTrain(refStart:length(spikeTrain)) = 0; 
    %        else 
    %            spikeTrain(refStart:refEnd) = 0; 
    %        end 
    %    end 
    % end 
    
   
    threshold = 0; 
    finalData = filteredData;
end 



end 