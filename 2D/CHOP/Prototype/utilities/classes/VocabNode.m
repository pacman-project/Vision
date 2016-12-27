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
        children@int32 % Label ids of this node's children in the previous level.
        modalExperts@int32; % Imagined modal reconstruction in terms of layer 1 nodes.
        categoryArr@single % 1 x numberOfCategories array. Sums up to 1. 
                    % Determines how much this node has been found in each
                    % category.
        minActivationLog@single; % Log of minimum activition in the training set.
    end
end
