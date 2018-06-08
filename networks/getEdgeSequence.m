function edgeSeq = getEdgeSequence(adjM)
%GETEDGESEQUENCE Summary of this function goes here
%   Detailed explanation goes here
%%% INPUT 
    % adjM : numNode x numNode x time-sequence, 3D matrix
    
% goal is to get this which works with plotDNarc
%  contactSequence = nEdges x 3 array of contacts (i,j,t)

edgeCount = 1; 
for t = 1:size(adjM, 3)
    edgeList = adj2edge(adjM(:, :, t)); 
    for e = 1:size(edgeList, 1)
        if ~isempty(edgeList)
            edgeSeq(edgeCount, 1) = edgeList(e, 1); 
            edgeSeq(edgeCount, 2) = edgeList(e, 2); 
            edgeSeq(edgeCount, 3) = t; 
            edgeCount = edgeCount + 1;
        end 
    end 
end 



end

