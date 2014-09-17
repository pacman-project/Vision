%> Class Name: GraphNode
%>
%> Description: At each level, realizations of the vocabulary parts is
%> stored in a GraphNode array. 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.11.2013
classdef GraphNode
    properties
        labelId % Label id of the part the realization belongs to.
        imageId % Image id of the realization. It can be used as the index 
                % in trainingFileNames array to get the name of the image.
        position % Position of the part's center in terms of image pixels. 
                 % Top left is (1,1).
        children % The children of this node in the previous level of 
                 % mainGraph(object graphs). They link to realization ids,
                 % which are indices of GraphNodes in respective mainGraph
                 % level.
        parents  % The parents of this node in the next level of mainGraph 
                 % (object graphs). Links to realization ids.
        adjInfo % Different than VocabNode, this field encodes relations of 
                % this graph node with its surroundings. It has nothing to
                % do with the children or the parent. After learning
                % statistical relations, the spatial information is stored
                % in this field which connects this node to other nodes in
                % its receptive field (thus forming receptive graph).
                % Is of the form [v1 v2 edgeId isDirected; ...]
        leafNodes % The leaf node ids (level 1 realizations) this node 
                  % corresponds to. Particularly important in visualization
                  % and calculation of coverage statistics.
        sign % Sign of the node. Assigned 0 unless the corresponding image 
             % is in background class.   
        dlReduction; % Reduction of DL of object graph due to this 
                     % instance.
    end
end