%> Name: getHistEdges
%>
%> Description: Given the nodes with their indices and center coordinates,
%> this function forms the edge structure depending on the relative
%> positions' placement on a pre-defined histogram matrix.
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
%> @retval edges Edges belonging to nodes.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 20.01.2014
function [ edges ] = getHistEdges( nodes, options, currentLevel, ~)
    %% Parameters
    edges = [];
    scale = (1/options.scaling)^currentLevel;
    neighborhood = options.edgeRadius * scale;
    nodeCount = size(nodes,1);
    
    %% Go over the node list and get the edges.
    nodeIds = cell2mat(nodes(:,1));
    nodeCoords = cell2mat(nodes(:,2));
    imageIds = cell2mat(nodes(:,3));
    
    for nodeItr = 1:nodeCount
        adjacentNodes = getadjacentNodes(nodeItr, nodes, neighborhood); 
        if numel(adjacentNodes)>0
            currEdges = [ones(size(adjacentNodes,1),1)*nodeItr, adjacentNodes, ones(size(adjacentNodes,1),1), zeros(size(adjacentNodes,1),1)];
            edges = [edges; currEdges];
        end
    end
end