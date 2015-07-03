%> Class Name: Substructure
%>
%> Description: The data structure in Subdue that encodes vocabulary nodes
%> (parts). Irrelevant fields for search are stripped out.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.11.2013
classdef Substructure
    properties
        centerId@int32 % The label id of the center node. For fast access.
        edges@int32 % Connectivity information of the children in this part's 
              % graph description.
        mdlScore@double % The MDL score of this node.
        normMdlScore@double % The normalized MDL score of this node. 
                     % (mdlScore/graphSize).
        instanceCenterIdx@int32 % An array of instance centers (realizations of the part)
        instanceChildren@int32 % An array of instance children nodes.
        instanceEdges@uint8; % The NxE array of instance edges, where N is the number of 
                                                   % realizations, and E is
                                                   % the number of edges in
                                                   % each insstance. Each number
                                                   % is an index to the
                                                   % array containing
                                                   % outgoing edges for
                                                   % each node.
       instanceSigns@uint8; % Nx1 sign array, 1 for foreground, 0 for background.
       instanceCategories@uint8; % Nx1 category array, category label for each instance.
       instanceMatchCosts@single; % Nx1 array, The matching cost of each instance to the sub.
       instanceValidationIdx@uint8; % Nx1 array, which tells which validation set 
                                    % this instance is drawn from.
    end
end