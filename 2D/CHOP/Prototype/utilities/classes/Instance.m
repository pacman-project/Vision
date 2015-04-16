%> Class Name: Instance
%>
%> Description: The data structure in Subdue that encodes mainGraph nodes
%> (realzations). Irrelevant fields for search are stripped out.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.11.2013
classdef Instance
    properties
        centerIdx@int32; % Center index (not label id!) of the center node 
                        % in this realization. The graph structure encoded
                        % in an Instance is isomorphic to the graph
                        % description of the Substructure this Instance
                        % belongs to. 
        edges@uint8; % Edges connecting the center node to the other nodes 
                    % in its receptive graph. Please note that this array
                    % does not encode all connections in the center node's
                    % receptive graph. It only retains the ones that are
                    % required to match this realization's graph
                    % instance to its linked part's graph description. 
        sign@uint8; % Signature of the instance, which typically should set 
                   % to 1 unless the relevant image is of background class.
        category@uint8; % Category of the the image this instance belongs to. ( This means max 255 categories are supported for now).
        matchCost@single; % Matching cost of this instance to the vocabulary part.
    end
end

