function spikeTrain = gmmDetect(data, windowLength, stepSize) 

    lowpass = 600; 
    highpass = 8000; 
    fs = 25000; 
    wn = [lowpass highpass] / (fs / 2); 
    filterOrder = 3;
    [b, a] = butter(filterOrder, wn); 
    filteredData = filtfilt(b, a, double(data));
    
   %  X = reshape(filteredData, windowLength, [])'; % reshape row-wise
    X = [ ]; % TODO: preallocate
    for cc = 1:windowLength 
        % X(:, cc) = filteredData(cc:(length(filteredData)-windowLength + cc - 1));
        column = round(linspace(1, length(filteredData) - windowLength +1, length(filteredData)/stepSize));
        X(:, cc) = filteredData(column + cc - 1); % increase stepsize to 10   
    end 
    
    options = statset('Display','final');
    gm = fitgmdist(X,2,'Options',options);
    % p = posterior(gm,X);
    idx = cluster(gm,X); % hard cluster
   
    
 
    

    %% Refractory period 

    % waveClus: 1.5ms 
    % Manuel's RCT paper: 2.0ms
    
    % shortcut, don't use in final implementation 
    if length(find(idx == 1)) < length(find(idx == 2))
        spikeTrain = double((idx == 1)); % check component is correct 
    else 
        spikeTrain = double((idx == 2)); 
    end 
   
    
    % fs = 25000;
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
    
end 