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
function [ nodes, parseTrees, nodeIds, nodeCoords ] = projectNode( nodes, vocabulary, ~, samplingMethod )
    levelItr = nodes(1,4);
    nodes = single(nodes);
    posDim = 2;
    parseTrees = zeros(size(nodes,1), levelItr);
    parseTrees(:,1) = (1:size(nodes,1))';
    nodeIds = zeros(size(nodes,1), levelItr);
    nodeIds(:,1) = nodes(:,1);
    nodeCoords = nodes(:,2:3);
    topLevel = levelItr;
    
    %% First, we recursively backproject the nodes. 
    nodeOffset = size(nodes,1);
    while levelItr > 1.001
        vocabLevel = vocabulary{levelItr};
        newNodes = cell(size(nodes,1),1);
        newParseTrees = cell(size(nodes,1),1);
        newNodeIds = cell(size(nodes,1),1);
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
                 for childItr = 2:numel(nodeCombination)
                    newNodeSet(childItr,2:3) = posVect(:,((childItr-2) * posDim + 1):((childItr-1)*posDim));
                 end
                 
                 % Shift nodes by an offset (mid of RF).
                 mins = min(newNodeSet(:, 2:3), [], 1);
                 maxs = max(newNodeSet(:,2:3), [], 1);
                 midPoint = round((mins+maxs)/2);
                 newNodeSet(:,2:3) = round(newNodeSet(:,2:3) - repmat(midPoint, size(newNodeSet,1),1));
            end
            % Finally, we update the positions by adding previous
            % offset.
            newNodeSet(:,2:3) = newNodeSet(:,2:3) + repmat(nodes(nodeItr,2:3), size(newNodeSet,1),1);
            
            % Assign rest of the fields and move on.
            newNodeSet(:, 1) = single(nodeCombination');
            newNodeSet(:, 4) = levelItr-1;
            newNodes{nodeItr} = newNodeSet;
            
            % Generate parse trees.
            tempParseTree = repmat(parseTrees(nodeItr,:), size(newNodeSet,1), 1);
            tempNodeIds = repmat(nodeIds(nodeItr,:), size(newNodeSet,1), 1);
            tempParseTree(:, (topLevel - levelItr) + 2) = ((nodeOffset+1):(nodeOffset+size(newNodeSet,1)))';
            tempNodeIds(:, (topLevel - levelItr) + 2) = newNodeSet(:,1);
            nodeCoords = [nodeCoords; newNodeSet(:,2:3)]; %#ok<AGROW>
            nodeOffset = nodeOffset + size(newNodeSet,1);
            newParseTrees{nodeItr} = tempParseTree;
            newNodeIds{nodeItr} = tempNodeIds;
        end
        
        % Shift the nodes by an offset.
        
        nodes = cat(1, newNodes{:});
        parseTrees = cat(1, newParseTrees{:});
        nodeIds = cat(1, newNodeIds{:});
        levelItr = levelItr - 1;
    end
    nodes = int32(round(nodes(:,1:3)));
    
    %% Then, we perform a simple inhibition process to remove overlapping level 1 instances.
    validNodes = ones(size(nodes,1),1) > 0;
    [~, firstSeenIdx] = unique(parseTrees, 'first');
    nodeIds = nodeIds(firstSeenIdx)';
    nodes = nodes(validNodes, :);
end

