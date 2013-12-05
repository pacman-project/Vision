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
%> @param options Program options.
%> @param currentLevel The level for which graph is extracted. Needed since
%> the local neighborhood is affected by level.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 20.11.2013
%> Ver 1.1 on 02.12.2013 'bin' property (edge type) added.
function [ edges ] = getEdges( nodes, options, currentLevel)
    %% Parameters
    edges = [];
    scale = (1/options.scaling)^currentLevel;
    neighborhood = options.edgeRadius * scale;
    
    nodeCount = size(nodes,1);
    
    for nodeItr = 1:nodeCount
        if strcmp(options.property, 'co-occurence')
            neighbors = getNeighbors(nodeItr, nodes, neighborhood); 
        end
        if numel(neighbors)>0
            currEdges = [ones(size(neighbors,1),1)*nodeItr, neighbors, ones(size(neighbors,1),1)];
            edges = [edges; currEdges];
        end
    end
    
end

%% Helper function to get neighbors of a node.
function [ neighbors ] = getNeighbors(nodeItr, nodes, neighborhood)
    centerArr = [ones(size(nodes,1), 1) * nodes{nodeItr,2}(1), ...
        ones(size(nodes,1), 1) * nodes{nodeItr,2}(2)];
    
    distances = sqrt(sum((centerArr - cell2mat(nodes(:,2))).^2, 2));
    neighbors = find(distances <= neighborhood);
    neighbors = neighbors(neighbors > nodeItr);
end
