%> Name: getEdges
%>
%> Description: Given the nodes with their indices and center coordinates,
%> this function forms the edge structure depending on the geometric property
%> to be examined.
%>
%> @param nodes The nodes extracted from the image in the following format:
%>      [ index x y;
%>        index x y;
%>          ...
%>        index x y]
%>              
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 20.11.2013
function [ edges ] = getEdges( nodes, property )
    %% Parameters
    edges = [];
    neighborhood = 5;
    
    nodeCount = size(nodes,1);
    
    for nodeItr = 1:nodeCount
        if strcmp(property, 'co-occurence')
            neighbors = getNeighbors(nodeItr, nodes, neighborhood); 
        end
        currEdges = [ones(size(neighbors,1),1)*nodeItr, neighbors];
        edges = [edges; currEdges];
    end
    
end

%% Helper function to get neighbors of a node.
function [ neighbors ] = getNeighbors(nodeItr, nodes, neighborhood)
    centerArr = [ones(size(nodes,1), 1) * nodes(nodeItr,2), ...
        ones(size(nodes,1), 1) * nodes(nodeItr,3)];
    
    distances = sqrt(sum((centerArr - nodes(:,[2 3])).^2, 2));
    neighbors = find(distances <= neighborhood);
    neighbors = neighbors(neighbors > nodeItr);
end
