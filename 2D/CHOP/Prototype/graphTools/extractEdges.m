%> Name: extractEdges
%>
%> Description: Given the node list, this function wraps the edge
%> extraction functions and returns the edges linking nodes calculated with
%> the correct method depending on options.property. There are two types of
%> edges: Object-wise edges consist of geometric relations between nodes
%> within a single object graph. Inter-object edges correspond to
%> high-level relations between realizations in different objects' graphs. 
%>
%> @param nodes The node list including label ids, positions and image ids
%> of each node.
%> @param mainGraph The object graphs' data structure.
%> @param options Program options.
%> @param currentLevel The currnet scene graph level.
%> @param datasetName Name of the dataset.
%> @param modes Modes up to level (currentLevel-1). If currentLevel == 0, 
%>      modes should be empty.
%> @param highLevelModes High level modes up to level (currentLevel-1). If 
%>      currentLevel == 0, similarly with modes, it should be empty.
%> 
%>
%> @retval modes The mode list representing edge categories.
%> @retval highLevelModes The high-level mode list representing inter-object-graph
%>      edge categories.
%> @retval edges Edges are of the form: [ node1, node2, label, directed;
%>                                        node1, node2, label, directed;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.12.2013
%> Addition of inter-image edges on 28.01.2014
function [modes, highLevelModes, edges, leafNodeAdjArr] = extractEdges(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName, modes, highLevelModes)
    %% Step 1: Learn low-level (within object) and high-level (between objects) modes.
    [currentModes, currentHighLevelModes] = learnStats(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName);
    if ~isempty(currentModes)
        modes = [modes, {currentModes}];
    end
    if ~isempty(currentHighLevelModes)
        highLevelModes = [highLevelModes, {currentHighLevelModes}];
    end

    %% Create within-object-graph edges.
    [edges, leafNodeAdjArr] = createEdgesWithLabels(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, modes);
    edges = unique(edges, 'rows', 'stable');
    
    %% Here, we create inter-object-graph edges. Nodes belonging to different object graphs are linked here.
    % CAUTION: Please note that no nodes within the SAME object graph should be
    % linked in this function. It ruins the inference process and causes it
    % to be unstable and inefficient.
    [highLevelEdges] = createHighLevelEdgesWithLabels(nodes, mainGraph, options, currentLevel, highLevelModes);
    highLevelEdges = unique(highLevelEdges, 'rows', 'stable');
    
    %% Combine both types of edges and return them.
    edges = [edges; highLevelEdges];
end

