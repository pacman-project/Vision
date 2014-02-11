%> Name: createHighLevelEdgesWithLabels
%>
%> Description: This function essentially analyzes relations between nodes
%> (realizations) in different object graphs. Please make sure that as a
%> result of this function, no edges between nodes within the same object's
%> graph are generated. Failure to do so may result in substantial
%> performance degredance and even hinder inference. 
%>
%> @param nodes The node list including label ids, positions and image ids
%>      of each node.
%> @param mainGraph The object graphs' data structure.
%> @param options Program options.
%> @param currentLevelId The current scene graph level id (1 .. n).
%> @param highLevelModes Cell array including high level modes for each
%>      level. Current level's modes can be accessed as
%>      highLevelModes{currentLevel}.
%> 
%> @retval edges Edges are of the form: [ nodeId1, nodeId2, edgeLabel, directed;
%>                                        nodeId3, nodeId4, edgeLabel, directed;
%>                                      ...]
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 28.01.2013, dummy function generation with proper
%>      input/output specs.
function [mainGraph] = createHighLevelEdgesWithLabels(mainGraph, options, currentLevel, highLevelModes)
    %% TODO : Implement high level edges here.
    if strcmp(options.highLevelProperty, 'multiview')
        % TODO: Implement 'multiview' type relations.
    elseif strcmp(options.highLevelProperty, 'category')
        % TODO: Implement 'category' type relations.
    else
        % 'none' option or invalid mode. No analysis on high level
        % relations.
    end
end
