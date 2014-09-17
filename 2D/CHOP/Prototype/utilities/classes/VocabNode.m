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
        label % Label id of the part. Same as the order of the part in the 
              % vocabulary array.
              % In case of redundant parts, the true label id of the part, 
              % indexed to vocabLevel, not redundantVocabLevel.
        mdlScore % Estimated MDL score of the part.
        normMdlScore % Normalized MDL score of the part.
        children % Label ids of this node's children in the previous level.
        parents  % Label ids of this node's parents in the next level.
        adjInfo  % Edge information of the form (v2 edgeId isDirected; 
                                               % v3 edgeId isDirected; ...)
        categoryArr % 1 x numberOfCategories array. Sums up to 1. 
                    % Determines how much this node has been found in each
                    % category.
    end
end