function gridTrace(electrodeMatrix, downFactor, removeElectrode, option, scaleb)
    % Practically Complete: 201801207 TS
    % plot the spike trace of the 60 electrodes 
    % can also take in 59 and 58, 
    % exect electrode matrix to be of dimension samples x nChannels
    % specify 8 x 8 grid 
    % TODO: use same axes scale (can be done with linkaxes)
    
    % Last update: TS 20180607
    % Author: Tim Sit, sitpakhang@gmail.com
    
    % INPUT 
    
    % downFactor 
        % the factor to downsample the traces by 
        % eg. downFactor = 1000, convert 25Khz to 25Hz 
        % if downFactor = 1, then no downsampling performed
    
    % removeElectrode
        % default = []
        % vector to remove electrode with specified numbers
        % note that you also have to include the grounded electrode if you
        % choose this option 
     
        
     % option  
        % default = [] 
        % if option is 'tight', then this uses tight subplot rather than
        % the matlab default
        
     % scaleb
        % default = 0
        % if scalebar is 1, then creates scalebar in the last subplot
        % note that this is currently not a very generisible option, ie. it
        % requires direct manipulation of the paramters within this
        % function.
        
    if ~exist('removeElectrode') 
        removeElectrode = [];
    end 
    
    if ~exist('option')
        option = []; 
    end 
    
     if ~exist('scaleb')
        scaleb = 0; 
    end 
    
        
    numRow = 8; 
    numColumn = 8; 
    trace = downsample(electrodeMatrix, downFactor);
    % pL = 1:(size(electrodeMatrix, 2)+3); % this plots in row-wise order
    pL = reshape(1:64, 8, 8)'; % this plots in column-wise order
    
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
        % pL = pL(pL ~= 5); 
        % pL = pL(pL ~= 16);
        pL = pL(pL ~= 33); 
        pL = pL(pL ~= 58); 
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
   
   %% Tight subplot 
   
   % utilitses tight subplot function 
   % avilable here: 
   % https://uk.mathworks.com/matlabcentral/fileexchange/27991-tight-subplot-nh--nw--gap--marg-h--marg-w-
   
   if strcmp(option, 'tight') 
       figure
       gap = [0.01, 0.01]; 
       marg_h = [0.05, 0.05]; 
       marg_w = [0.01, 0.01]; 
       ha = tight_subplot(numRow, numColumn, gap, marg_h, marg_w);
       
       % remove the corners created in ha 
       axes(ha(1)); removeAxis; 
       axes(ha(8)); removeAxis;
       axes(ha(57)); removeAxis;
       axes(ha(64)); removeAxis;
       
       for plotN = 1:size(electrodeMatrix, 2)
           axes(ha(pL(plotN))); 
           plot(trace(:, plotN))
           % title(num2str(plotN))
           xlim([1 length(trace)]) % this makes the plot xlim proper!
           % title(num2str(plotN))
           aesthetics()
           removeAxis()
           lineThickness(3)
           hold on
       end 
       % linkaxes() % to make them have the same x and y limits
       % use the grounded electrode (number 4), to set y limits
       % ylim([min(electrodeMatrix(:, 4)), max(electrodeMatrix(:, 4))])
       % ylim([min(min(electrodeMatrix)), max(max(electrodeMatrix))])
       % make it square-ish
       set(gcf, 'Position', [10, 10, 1500, 800])
       linkaxes()
   end 
   
   %% Scalebar 
   
   if scaleb == 1 && ~strcmp(option, 'tight')
       subplot(numRow, numColumn, 64)
       linkaxes
       sb = scalebar;
       sb.YLen = 50; 
       sb.XLen = 12500; 
       sb.YUnit = '\muV';
       sb.XUnit = 'ms'; 
       sb.Position = [6000, -80];
   elseif scaleb == 1 && strcmp(option, 'tight')
       axes(ha(64))
       sb = scalebar;
       % sb.YLen = 50; 
       % sb.XLen = 12500; 
       sb.YUnit = '\muV';
       sb.XUnit = 'ms'; 
       % sb.Position = [8000, -60];
       sb.YLen = 30; sb.XLen = 50; sb.Position = [110 -20]
   end 
   
end 