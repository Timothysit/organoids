function plotAdj(adjM, goodElectrodes)
% note that this is quite customised for the mecp2 project (MEA 8 x 8) data
% INPUT 
    % plotAdj | N x N adjacency matrix (where N is the number of active electrodes) 
    % goodElectrodes | N x 1 vector saying which electrodes are active;
    % only plot those

% Network adjacency plot, like the one by Manuel 
% author: Tim Sit github/timothysit 
% last update: 20180407 
% TODO: check that the coordinates agrees with other heatmap visualisations

%% generate the possible 2-val combinations of 1 to 8, with repetition

% n = 8; k = 2;
% nk = nmultichoosek(1:n,k);
% coord=zeros(0,k);
% for i=1:size(nk,1)
%     pi = perms(nk(i,:));
%     coord = unique([coord; pi],'rows');
% end
% 
% coord(:, 1) = 9 - coord(:, 1); 

% coord needs to count from top to bottom columnwise, ie. [8, 1], [8,
% 2]... 

yTemp = 1:8; 
yCoord = repmat(fliplr(yTemp), 1, 8);

xTemp = 1:8; 
xCoord = repelem(xTemp, 8);

coord = [xCoord', yCoord'];

%% prune coord 


coord(coord(:, 1) == 1 & coord(:, 2) == 1, :) = [];
coord(coord(:, 1) == 1 & coord(:, 2) == 8, :) = [];
coord(coord(:, 1) == 8 & coord(:, 2) == 8, :) = [];
coord(coord(:, 1) == 8 & coord(:, 2) == 1, :) = [];

coord = coord(goodElectrodes, :);





%% do some pruning of the adjacency so we don't plot everytying 

edges = adjM; % flip things to make it agree with the controllabiliyt heatmap
% edges(edges == 1) = 0; % remove self-correlation 
% threshold = prctile(edges(:), 75); % get the 95th percentile 
threshold = 0;
edges(edges < threshold) = 0.0001; 
edges(isnan(edges)) = 0.0001; % I think 0 doens't quite work


%% actual network plot 
% figure
% gplot(a,coord) % most simple way of doing it
% h = adjacency_plot_und(a, coord); % requires brain connectivity toolbox

weights = sum(edges) / sum(sum(edges));
wgPlot(edges, coord, 'edgeColorMap', colormap(flipud(gray(256))), 'vertexWeight', weights', 'vertexScale',200,'edgeWidth',2);
% caxis([0 1])
% vertex scale oringal value: 200
% edgewidth original value: 2

aesthetics

end 