%> Name: projectNode
%>
%> Description: Given the nodes in 'nodes' array, and vocabulary descriptions 
%> in vocabulary, we backproject the nodes to image plane by recursively
%> projecting sub-nodes.
%>
%> @param nodes exported nodes in the format of 
%> [labelId, centerPosX, centerPosY, levelId]
%>   
%> @param vocabulary The vocabulary.
%> @param inhibitionRadius The radius in which we will suppress other responses.
%>  
%> @retval nodes Projected level 1 nodes. 
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 12.08.2015
function [ nodes, subChildrenExperts, subChildren, orNodeChoices, orNodeChoiceCounts ] = projectNode( nodes, vocabularyDistributions, samplingMethod, options, orNodeChoice )

    if nargin == 4
         orNodeChoice = [];
    end
    
    minOptimizationLayer = 3;
    levelItr = nodes(1,4);
    nodes = single(nodes);
    posDim = 2;
    subChildren = cell(size(nodes,1), 1);
    subChildrenExperts = cell(size(nodes,1), 1);
    orNodeChoices = zeros(size(nodes,1), 1);
    orNodeChoiceCounts = zeros(size(nodes,1), 1);
    topLevel = levelItr;
    
    %% First, we recursively backproject the nodes. 
    startNodeOffset = 0;
    while levelItr > 1.001
        % Get receptive field size.
        rfSize = getRFSize(options, levelItr);
        centerRF = floor(rfSize(1)/2)+1;
         
        % 
        if levelItr < minOptimizationLayer
             samplingMethod = 'modal';
        end
        vocabLevelDistributions = vocabularyDistributions{levelItr};
        prevLevelDistributions = vocabularyDistributions{levelItr-1};
        newNodes = cell(size(nodes,1),1);
        
        % If previous layer already has modal experts defined, we don't
        % need to proceed.
        stopFlag = ~isempty(prevLevelDistributions(1).modalExperts) && strcmp(samplingMethod, 'modal') && levelItr == topLevel;
        
        for nodeItr = 1:size(nodes,1)
            vocabNodeDistributions = vocabLevelDistributions(nodes(nodeItr,1));
            newNodeSet = zeros((size(vocabNodeDistributions.childrenLabelDistributions,2)-1), 4, 'single');
            
            % Sample from the discrete and continuous distributions.
            childrenLabelDistributions = vocabNodeDistributions.childrenLabelDistributions;
            choiceCounts = size(childrenLabelDistributions,1);

            % Reconstruct the nodes.
            if ~isempty(orNodeChoice)
                 assignedRow = orNodeChoice;
            elseif strcmp(samplingMethod, 'modal')
                 [~, assignedRow] = max(childrenLabelDistributions(:, end));
            else
                 % Sample from the discrete distribution. We obtain a
                 % random number and assign it to an interval to calculate
                 % probabilities. 
                 probabilities = childrenLabelDistributions(:,end);
                 probabilities = [0; probabilities]; %#ok<AGROW>
                 cumProbabilities = cumsum(probabilities);
                 cumProbabilities(end) = cumProbabilities(end) + 0.0001;
                 cumProbabilities(1) = cumProbabilities(1) - 0.0001;
                 randNumber = single(rand());
                 
                 % Find which row to assign.
                 assignedRow = find(randNumber >= cumProbabilities(1:(end-1))  & randNumber < cumProbabilities(2:end));
            end
            nodeCombination = childrenLabelDistributions(assignedRow, 1:(end-1)); 
            
            % Given the node combinations, we obtain relevant position
            % distributions and sample from that. Mode of a Gaussian
            % mixture is obtained by simply getting the distribution with
            % most contribution, and sampling from that. 
            % TODO: We are planning to replace this with a smarter
            % search mechanism.
            posDistributions = vocabNodeDistributions.childrenPosDistributions;
            posDistributionModes = vocabNodeDistributions.childrenPosDistributionModes;
 
            if ~isempty(posDistributions.mu)
                 % Sample from the distribution.
                 if strcmp(samplingMethod, 'modal')
                     assignedDist = posDistributionModes(assignedRow);
                     posVect = posDistributions.mu(assignedDist,:);
                 else
                     posVect = random(posDistributions, 1);
                 end
                 
                 % Assign positions.
                 for childItr = 1:numel(nodeCombination)
                    newNodeSet(childItr,2:3) = posVect(:,((childItr-1) * posDim + 1):((childItr)*posDim));
                 end
%                  
%                  % Shift nodes by an offset (mid of RF).
%                  mins = min(newNodeSet(:, 2:3), [], 1);
%                  maxs = max(newNodeSet(:,2:3), [], 1);
%                  midPoint = round((mins+maxs)/2);
%                  newNodeSet(:,2:3) = round(newNodeSet(:,2:3) - repmat(midPoint, size(newNodeSet,1),1));
            end
            
            % Finally, we update the positions by adding previous
            % offset.
            newNodeSet(:,2:3) = newNodeSet(:,2:3) + repmat(nodes(nodeItr,2:3), size(newNodeSet,1),1);
            newNodeSet(:, 1) = single(nodeCombination');
            newNodeSet(:, 4) = levelItr-1;
            
            % Assign rest of the fields and move on.
            newNodes{nodeItr} = newNodeSet;
            
             % If there are existing predictions in previous layer, let's
             % combine them.
             if levelItr > 1 && stopFlag
                    % Combine experts with the imaginations from the previous
                    % layer.
                    newNodeSet = int32(round(newNodeSet(:,1:3)));
                    allExperts = cell(size(newNodeSet,1),1);
                    for expertItr = 1:size(newNodeSet,1)
                         newChildren = prevLevelDistributions(newNodeSet(expertItr,1)).modalExperts;
                         newChildren(:,2:3) = repmat(newNodeSet(expertItr,2:3), size(newChildren,1),1) + newChildren(:,2:3);
                         allExperts{expertItr} = newChildren;
                    end
                    if levelItr == topLevel
                         subChildrenExperts{nodeItr} = allExperts;
                    end
             end
            
            % Save this choice.
             if topLevel == levelItr
                  subChildren(nodeItr) = {newNodeSet};
                  orNodeChoices(nodeItr) = assignedRow;
                  orNodeChoiceCounts(nodeItr) = choiceCounts;
             end
        end
        
        % If stop flag is on, we don't need to proceed further.
        if stopFlag
             nodes = cat(1, subChildrenExperts{:});
             nodes = cat(1, nodes{:});
             break;
        end
        
        % Shift the nodes by an offset.
        startNodeOffset = size(nodes,1) + startNodeOffset;
        nodes = cat(1, newNodes{:});
        levelItr = levelItr - 1;
    end
    
    nodes = int32(round(nodes(:,1:3)));
end

