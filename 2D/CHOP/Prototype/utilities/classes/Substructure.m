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
        instances % An array of instances (realizations of the part).
    end
end