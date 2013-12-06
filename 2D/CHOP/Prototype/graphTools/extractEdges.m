%> Name: extractEdges
%>
%> Description: Given the node list, this function wraps the edge
%> extraction functions and returns the edges linking nodes calculated with
%> the correct method depending on options.property.
%>
%> @param nodes The node list including label ids, positions and image ids
%> of each node.
%> @param options Program options.
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
function [modes, edges] = extractEdges(nodes, options, currentLevel)
    modes = [];
    if strcmp(options.property, 'mode')
        % If needed, process nodes to determine 'mode's, i.e. pair-wise relations.
        [modes, edges] = addModes(nodes, options, currentLevel);
    else
        % Form simple edges (co-occurence based)
        edges = getEdges(nodes, options, currentLevel);
    end
end

