%> Class Name: Substructure
%>
%> Description: The data structure in Subdue that encodes vocabulary nodes
%> (parts). Irrelevant fields for search are stripped out. This class has
%> the instances of subs in it as well. This method proved to be less costly
%> than representing each instance as a single entity, and keeping an array
%> of instances. 
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
       instanceCenterIdx@int32 % An array of instance centers (realizations of the part)
       instanceChildren@int32 % An array of instance children nodes.
       instanceSigns@uint8; % Nx1 sign array, 1 for foreground, 0 for background.
    end
end