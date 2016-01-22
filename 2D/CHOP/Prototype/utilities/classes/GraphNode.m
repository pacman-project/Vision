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
        labelId@int32 % Label id of the part the realization belongs to.
        nodeProbability@single % The probability of the node label.
        realLabelId@int32 % Actual label of the associated part (considering OR nodes)
        imageId@int32 % Image id of the realization. It can be used as the index 
                % in trainingFileNames array to get the name of the image.
        position@int32 % Position of the part's center in terms of image pixels, 
                       % Keep in mind that the resolution is usually halved at each level. 
                 % Top left is (1,1).
        precisePosition@single % Exact of the part's center in terms of 
                 % continuous image coordinates. % Top left is (1,1).
        children@int32 % The children of this node in the previous level of 
                 % mainGraph(object graphs). They link to realization ids,
                 % which are indices of GraphNodes in respective mainGraph
                 % level.
        mapping@uint8 % The mapping of children nodes to the definition of the
                      % part, indicated by labelId. 
        parents@int32  % The parents of this node in the next level of mainGraph 
                 % (object graphs). Links to realization ids.
        adjInfo@int32 % Different than VocabNode, this field encodes relations of 
                % this graph node with its surroundings. It has nothing to
                % do with the children or the parent. After learning
                % statistical relations, the spatial information is stored
                % in this field which connects this node to other nodes in
                % its receptive field (thus forming receptive graph).
                % Is of the form [v1 v2 edgeId isDirected; ...]
        edgeProbabilities@single % The probability of each edge occuring given the mode.
        realEdgeLabels@int32 % The true labels of each edge.
        leafNodes@int32 % The leaf node ids (level 1 realizations) this node 
                  % corresponds to. Particularly important in visualization
                  % and calculation of coverage statistics.
        activation@single; % Used in inference. Signals our confidence in this node.
        likelihood@single; % log likelihood of this particular node. 
        sign % Sign of the node. Assigned 0 unless the corresponding image 
             % is in background class.   
    end
end