function plotDynamicAdj(spikes) 
% Plot dynamic adjacency matrix 
% Requires Dynamic-Graph-Metrics-master 

% start with spikes 
% end with plot

adjM = dynamicAdjacency(spikes, 50, 'covariance'); 
% edgeList = adj2edge(adjM); 
edgeSeq = getEdgeSequence(adjM); 
plotDNarc(edgeSeq, size(spikes, 2)); 

aesthetics(); 

end

