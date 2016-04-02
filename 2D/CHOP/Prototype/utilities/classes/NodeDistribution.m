%> Class Name: NodeDistribution
%>
%> Description: Only preserves relevant generation information. 
%> Used for part reconstruction.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 30.03.2016
classdef NodeDistribution
    properties
        childrenLabelDistributions@single % Array containing possible label combinations
                                                                 % and their probabilities.
        childrenPosDistributions@gmdistribution; % Cell array of multi-modal gaussian distributions 
                                                                 % modelling joint space of children positions. 
                                                                 % One for every discrete combination.
        childrenPosDistributionModes@uint8 % Vector that has as many elements as the number of rows in 
                                                                    % childrenLabelDistributions. Each element shows which mode of the
                                                                    % distribution in childrenPosDistributions is used for
                                                                    % reconstructing each label combination.
        modalExperts@int32; % Imagined modal reconstruction in terms of layer 1 nodes.
                                            % Used for fast processing.
    end
end
