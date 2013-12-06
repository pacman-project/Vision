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
function [ edges ] = getEdges( nodes, options, currentLevel)
    %% Parameters
    edges = [];
    scale = (1/options.scaling)^currentLevel;
    neighborhood = options.edgeRadius * scale;
    
    nodeCount = size(nodes,1);
    
    %% Go over the node list and get the edges.
    for nodeItr = 1:nodeCount
        adjacentNodes = getadjacentNodes(nodeItr, nodes, neighborhood); 
        if numel(adjacentNodes)>0
            currEdges = [ones(size(adjacentNodes,1),1)*nodeItr, adjacentNodes, ones(size(adjacentNodes,1),1), zeros(size(adjacentNodes,1),1)];
            edges = [edges; currEdges];
        end
    end
    
end

%% Helper function to get adjacentNodes of a node.
function [ adjacentNodes ] = getadjacentNodes(nodeItr, nodes, neighborhood)
    centerArr = [ones(size(nodes,1), 1) * nodes{nodeItr,2}(1), ...
        ones(size(nodes,1), 1) * nodes{nodeItr,2}(2)];
    
    distances = sqrt(sum((centerArr - cell2mat(nodes(:,2))).^2, 2));
    adjacentNodes = find(distances <= neighborhood);
    adjacentNodes = adjacentNodes(adjacentNodes > nodeItr);
end

%% Helper function to estimate relative positioning of two nodes.
function [vector] = getPositionVector(node1, node2, nodes, neighborhood)
    point1 = nodes{node1,2};
    point2 = nodes{node2,2};
    vector = (point2 - point1) / neighborhood;
end

