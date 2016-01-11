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
function [ nodes ] = projectNode( nodes, vocabulary, inhibitionRadius, samplingMethod )
    levelItr = nodes(1,4);
    nodes = single(nodes);
    inhibitionRadius = inhibitionRadius - 1;
    posDim = 2;
    
    %% First, we recursively backproject the nodes. 
    while levelItr > 1.001
        vocabLevel = vocabulary{levelItr};
        newNodes = cell(size(nodes,1),1);
        for nodeItr = 1:size(nodes,1)
            vocabNode = vocabLevel(nodes(nodeItr,1));
            newNodeSet = zeros(numel(vocabNode.realChildren), 4, 'single');
            
            % Sample from the discrete and continuous distributions.
            childrenLabelDistributions = vocabNode.childrenLabelDistributions;
            if strcmp(samplingMethod, 'modal')
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
            posDistributions = vocabNode.childrenPosDistributions{1};
            posDistributions = posDistributions{assignedRow};
            if ~isempty(posDistributions)
                 mixturePs = (posDistributions.PComponents)'; 
                 % Sample from the distribution.
                 if strcmp(samplingMethod, 'modal')
                     [~, assignedDist] = max(mixturePs); %#ok<UDIM>
                     posVect = posDistributions.mu(assignedDist,:);
                 else
                     posVect = random(posDistributions, 1);
                 end
                 
                 % Assign positions.
                 newNodeSet(1,2:3) = nodes(nodeItr,2:3);
                 for childItr = 2:numel(nodeCombination)
                    newNodeSet(childItr,2:3) = posVect(:,((childItr-2) * posDim + 1):((childItr-1)*posDim)) + nodes(nodeItr,2:3);
                 end
            end
            
            % Assign rest of the fields and move on.
            newNodeSet(:, 1) = single(nodeCombination');
            newNodeSet(:, 4) = levelItr-1;
            newNodes{nodeItr} = newNodeSet;
        end
        nodes = cat(1, newNodes{:});
        levelItr = levelItr - 1;
    end
    nodes = int32(round(nodes(:,1:3)));
    
    %% Then, we perform a simple inhibition process to remove overlapping level 1 instances.
    validNodes = ones(size(nodes,1),1) > 0;
    nodes = nodes(validNodes, :);
end

