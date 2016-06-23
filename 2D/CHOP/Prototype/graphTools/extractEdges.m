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
%>
%> @retval edges Edges are of the form: [ node1, node2, label, directed;
%>                                        node1, node2, label, directed;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.12.2013
%> Addition of inter-image edges on 28.01.2014
function [graphLevel] = extractEdges(graphLevel, firstLevelAdjNodes, options, currentLevelId)
    %% Create within-object-graph edges.
    edgeNoveltyThr = max(options.minEdgeNoveltyThr, options.edgeNoveltyThr - options.edgeNoveltyThrRate * max(0, (currentLevelId-2)));
    display('Extracting edges...');
    if ~isempty(graphLevel)
        [graphLevel] = createEdgesWithLabels(graphLevel, firstLevelAdjNodes, options, currentLevelId, edgeNoveltyThr);
    end
end