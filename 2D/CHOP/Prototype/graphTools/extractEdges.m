%> Name: extractEdges
%>
%> Description: Given the node list, this function wraps the edge
%> extraction functions and returns the edges linking nodes calculated with
%> the correct method depending on options.property. There are two types of
%> edges: Object-wise edges consist of geometric relations between nodes
%> within a single object graph. Inter-object edges correspond to
%> high-level relations between realizations in different objects' graphs. 
%>
%> @param mainGraph The object graphs' data structure.
%> @param options Program options.
%> @param currentLevelId The current scene graph level id.
%> @param modes Modes up to level (currentLevel-1). If currentLevel == 0, 
%>      modes should be empty.
%> @param highLevelModes High level modes up to level (currentLevel-1). If 
%>      currentLevel == 0, similarly with modes, it should be empty.
%> 
%>
%> @retval modes The mode list representing edge categories.
%> @retval edges Edges are of the form: [ node1, node2, label, directed;
%>                                        node1, node2, label, directed;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.12.2013
%> Addition of inter-image edges on 28.01.2014
function [modes, mainGraph] = extractEdges(mainGraph, options, currentLevelId, modes)
    %% Step 1: Learn low-level (within object) and high-level (between objects) modes.
    if options.isTraining
        [currentModes] = learnStats(mainGraph, options, currentLevelId);
        if ~isempty(currentModes)
            modes = [modes, {currentModes}];
        end
    end
    load('hMatrix.mat', 'hMatrix'); 
    %% Create within-object-graph edges.
    display('Extracting edges...');
    [mainGraph] = createEdgesWithLabels(mainGraph, options, currentLevelId, modes, hMatrix);
end

