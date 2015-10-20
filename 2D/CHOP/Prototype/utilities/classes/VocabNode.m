%> Class Name: VocabNode
%>
%> Description: At each level, the vocabulary parts are stored in a
%> VocabNode array. 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.11.2013
classdef VocabNode
    properties
        label@int32 % Label id of the part. Same as the order of the part in the 
              % vocabulary array.
              % In case of redundant parts, the true label id of the part, 
              % indexed to vocabLevel, not redundantVocabLevel.
        mdlScore@double % Estimated MDL score of the part.
        normMdlScore@double % Normalized MDL score of the part.
        children@int32 % Label ids of this node's children in the previous level.
        parents@int32 % Label ids of this node's parents in the next level.
        adjInfo@int32 % Edge information of the form (v2 edgeId isDirected; 
                                               % v3 edgeId isDirected; ...)
        childrenNodeProbs@single % Nx2 array, in the form of [realNodeId, probability; ...]
                                                     % which assigns a probability to each
                                                     % alternative or node label.
        childrenPosMean@single % The relative position of the children with 
                           % respect to the center of gravity. Used to
                           % perform precise reconstruction of nodes. (Nx2,
                           % with each row representing x, y position of a
                           % child.
        childrenPosCov@single % The relative position of the children with 
                           % respect to the center of gravity. Used to
                           % perform precise reconstruction of nodes. (Nx2,
                           % with each row representing x, y position of a
                           % child.
        categoryArr@single % 1 x numberOfCategories array. Sums up to 1. 
                    % Determines how much this node has been found in each
                    % category.
        numberOfProbabilisticChoices@int32; % Includes the number of 
                     % probabilistic choices we need to make to match this
                     % node to a chunk of data across the hierarchy.
        orgRank@int32;
    end
end
