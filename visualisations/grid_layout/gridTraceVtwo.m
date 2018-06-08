function gridTrace(electrodeMatrix, downFactor, removeElectrode)
    % Practically Complete: 20180127 TS
    % This V2 version is slightly modified for the organoid project
    % plot the spike trace of the 60 electrodes 
    % can also take in 59 and 58, 
    % exect electrode matrix to be of dimension samples x nChannels
    % specify 8 x 8 grid 
    % TODO: use same axes scale (can be done with linkaxes)
    
    % Last update: TS 20180426
    
    % INPUT 
    
        % removeElectrode, default = []
        % vector to remove electrode with specified numbers
        % note that you also have to include the grounded electrode if you
        % choose this option 
    if ~exist('removeElectrode') 
        removeElectrode = [];
    end 
    
    
        
    numRow = 8; 
    numColumn = 8; 
    trace = downsample(electrodeMatrix, downFactor);
    pL = 1:(size(electrodeMatrix, 2)+3); % this plots in row-wise order
    % pL = reshape(1:64, 8, 8)'; % this plots in column-wise order
    
    %% Absent electrodes
    
    pL = pL(pL~=1); 
    pL = pL(pL~=8); 
    pL = pL(pL~=57);
    
    %% Grounding Electrodes
    % note that part is specifically for the MEA Project 2017-2018
    if ~isempty('removeElectrode') 
        % do nothing
    elseif size(electrodeMatrix, 2) == 59 
        % pL = pL(pL ~= 5); % for horizontal, which is wrong
        pL = pL(pL ~= 33); % for vertical, which is right
    elseif size(electrodeMatrix, 2) == 58
        pL = pL(pL ~= 5); 
        pL = pL(pL ~= 16);
        % pL = pL(pL ~= 33); % this is for vertical (HPC)
        % pL = pL(pL ~= 58); % this is for vertical (HPC dataset)
    end 
    
    %% Custom removed electroeds 
    
    if ~isempty(removeElectrode) 
        for i = 1:length(removeElectrode)
            pL = pL(pL ~= removeElectrode(i) );
        end 
    end 
    
    
    
   %% Actual Plot 
   for plotN = 1:size(electrodeMatrix, 2)
       subplot(numRow, numColumn, pL(plotN))
       plot(trace(:, plotN))
       xlim([1 length(trace)]) % this makes the plot xlim proper!
       title(num2str(plotN))
       aesthetics()
       removeAxis()
       hold on
   end 
   % linkaxes() % to make them have the same x and y limits
   % use the grounded electrode (number 4), to set y limits
   % ylim([min(electrodeMatrix(:, 4)), max(electrodeMatrix(:, 4))])
   % ylim([min(min(electrodeMatrix)), max(max(electrodeMatrix))])
   % make it square-ish
   set(gcf, 'Position', [10, 10, 1500, 800])
   linkaxes()
   %suptitle('1209 DIV22')
   
   % scalebar (didn't manage to get to work properly)
   % sb = scalebar;
   % sb.Position = [55, -10];
   % sb.XLen = 15;  
   % sb.YLen = 2;
   
end 