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
        centerIdx = []; % Center index (not label id!) of the center node 
                        % in this realization. The graph structure encoded
                        % in an Instance is isomorphic to the graph
                        % description of the Substructure this Instance
                        % belongs to. 
        edges = []; % Edges connecting the center node to the other nodes 
                    % in its receptive graph. Please note that this array
                    % does not encode all connections in the center node's
                    % receptive graph. It only retains the ones that are
                    % required to match this realization's graph
                    % instance to its linked part's graph description. Is
                    % of the form (v2, edgeId, isDirected; 
                    % v3, edgeId, isDirected; ...]
        sign = []; % Signature of the instance, which typically should set 
                   % to 1 unless the relevant image is of background class.
        dlReduction = []; % Reduction of DL of object graph due to this 
                   % instance.
    end
end

