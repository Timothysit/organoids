function edgeList = adj2edge(adj)
%ADJ2EDGE convert binary adjacency matrix to edgeList
%   %%% INPUT 
          % adj : symmetrical, square, matrix with 1 and 0
          
[r, c] = find(adj); 
edges = [r, c]; 

nonSelfEdges = edges(find(r-c ~=0), :); 

edgeList = unique(nonSelfEdges, 'rows');

end

