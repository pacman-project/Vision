%> Name: extractEdges
%>
%> Description: Given the node list, this function wraps the edge
%> extraction functions and returns the edges linking nodes calculated with
%> the correct method depending on options.property.
%>
%> @param nodes The node list including label ids, positions and image ids
%> of each node.
%> @param mainGraph The object graphs' data structure.
%> @param options Program options.
%> @param currentLevel The currnet scene graph level.
%> @param datasetName Name of the dataset.
%> @param modes If empty, new modes are to be learned. If not, edge labels
%> will be formed depending on existing modes.
%>
%> @retval modes The mode list representing edge categories.
%> @retval edges Edges are of the form: [ node1, node2, mode, directed;
%>                                        node1, node2, mode, directed;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.12.2013
function [modes, edges, leafNodeAdjArr] = extractEdges(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName, modes)
    if strcmp(options.property, 'mode')
        % If needed, process nodes to determine 'mode's, i.e. pair-wise relations.
        [modes, edges, leafNodeAdjArr] = calculateAndAssignModes(nodes, mainGraph, leafNodeAdjArr, options, currentLevel, datasetName, modes);
    else
        % Form simple edges (co-occurence based)
        modes = [];
        leafNodeAdjArr = [];
        edges = getEdges(nodes, options, currentLevel, datasetName);
    end
end

