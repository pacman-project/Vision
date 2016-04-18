classdef TreeNode
     %UNTITLED Summary of this class goes here
     %   Detailed explanation goes here
     
     properties
           realLabelId@int32
           precisePosition@int32
           subChildrenExperts@int32
           subChildren@TreeNode
           orNodeChoices@double
           orNodeChoiceCounts@double
     end
end