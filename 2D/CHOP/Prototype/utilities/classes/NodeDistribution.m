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
        minPosActivationLog@single; % Min position activation.
        childrenPosDistributionModes@uint8 % Vector that has as many elements as the number of rows in 
                                                                    % childrenLabelDistributions. Each element shows which mode of the
                                                                    % distribution in childrenPosDistributions is used for
                                                                    % reconstructing each label combination.
        modalExperts@int32; % Imagined modal reconstruction in terms of layer 1 nodes.
                                            % Used for fast processing.
    end
    
    methods
         function [partialLabels, instanceJointPos] = predictMissingInfo(obj, partialLabels, partialChildren, precisePositions)
             numberOfInstances=  size(partialLabels,1);
             numberOfChildren = size(partialLabels,2);
             muArr = obj.childrenPosDistributions.mu;
             instanceJointPos = zeros(numberOfInstances, (numberOfChildren-1)*2, 'single');
             
             % Go over instances and 
             for instanceItr = 1:numberOfInstances
                  relevantLabels = partialLabels(instanceItr,:);
                  centerCoords = precisePositions(partialChildren(instanceItr,1),:);
                  
                  % Obtain closest label combination.
                  validCols = [~isinf(relevantLabels), true];
                  invalidChildCols = ~validCols(1:end-1);
                  reducedChildrenLabelDistributions = obj.childrenLabelDistributions(:, validCols);
                  [~, idx] = ismember(relevantLabels(validCols(1:end-1)), reducedChildrenLabelDistributions(:,1:end-1), 'rows');
                  
                  % Based on the partial data, let's pick up the label
                  % combinations.
                  if isempty(idx) || idx == 0
                       probs = obj.childrenLabelDistributions(:, end);
                       [~, assignedRow] = max(probs);
                  elseif numel(idx) > 1
                       assignedRow = idx;
                  else
                       probs = obj.childrenLabelDistributions(idx, end);
                       [~, idx2] = max(probs);
                       assignedRow = idx(idx2);
                  end
                  
                  % Complete missing node labels here.
                  partialLabels(instanceItr, invalidChildCols) = obj.childrenLabelDistributions(assignedRow, ~validCols);
                  
                  % Now, we find closest locations for the nodes.
                  jointPos = zeros(1, numberOfChildren * 2);
                  secChildren = partialChildren(instanceItr, :);
                  for childItr = 1:numel(secChildren)
                       % Missing children will be Inf!
                       if isinf(secChildren(childItr))
                            jointPos(:, (childItr-1)*2+1:(childItr)*2) = Inf;
                       else
                            jointPos(:, (childItr-1)*2+1:(childItr)*2) = precisePositions(secChildren(childItr), :) - centerCoords;
                       end
                  end
                  instanceJointPos(instanceItr, :) = jointPos;
                  
                  % We have calculated joint position. Now, we'll obtain
                  % closest mode.
                  invalidCols = isinf(jointPos);
                  validCols = ~invalidCols;
                  if nnz(validCols) < numel(validCols)
                       if nnz(validCols) == 0
                            assignedMode = obj.childrenPosDistributionModes(assignedRow);
                       elseif nnz(validCols) ~= numel(validCols)
                            % Finally, calculate the final position.
                            shortMuArr = muArr(:, validCols);
                            [~, assignedMode] = pdist2(shortMuArr, jointPos(:, validCols), 'euclidean', 'Smallest', 1);
                       end
                       instanceJointPos(instanceItr, invalidCols) = round(muArr(assignedMode, invalidCols));
                  end
             end
         end
    end
end