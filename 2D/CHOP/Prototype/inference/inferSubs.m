%> Name: runSubdue
%>
%> Description: Inference on graphLevel with given subs in vocabLevel. Each
%> composition is searched for in graphLevel.
%>
%> @param vocabLevel Input vocabulary level. Compositions in this vocabulary 
%> level are detected in graphLevel.
%> @param graphLevel The current object graphs' level. The graphLevel's
%> nodes are sorted first by their imageId, then labelId.
%> @param options Program options.
%>
%> @retval graphLevel The graph level consisting of discovered nodes.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 05.02.2014
function [exportArr, activationArr, allPrecisePositions] = inferSubs(vocabulary, nodes, allModes, nodeActivations, distanceMatrices, optimalThresholds, edgeChangeLevel, options)
    % Read data into helper data structures.
    edgeDistanceMatrix = double(options.edgeDistanceMatrix);
    firstLevelAdjNodes = [];
    
    % If fast inference is not required, we do not perform inhibition.
    singlePrecision = options.singlePrecision;
    
    exportArr = [];
    if isempty(nodes) || isempty(vocabulary)
        return;
    end
    allNodes = cell(numel(vocabulary),1);
    allPrecisePositions = cell(numel(vocabulary), 1);
    allActivations = cell(numel(vocabulary),1);
    allNodes(1) = {nodes};
    allActivations(1) = {nodeActivations};
    leafNodeArr = num2cell((int32(1:size(nodes,1)))');
    precisePositions = single(nodes(:,4:5));
    allPrecisePositions(1) = {precisePositions};
    
    for vocabLevelItr = 2:numel(vocabulary)
        poolDim = options.poolDim;
        
        %% Match subs from vocabLevel to their instance in graphLevel.
        vocabLevel = vocabulary{vocabLevelItr};
        vocabRealizations = cell(numel(vocabLevel),1);
        vocabRealizationsActivations = cell(numel(vocabLevel),1);
        nodeDistanceMatrix = distanceMatrices{vocabLevelItr-1};
        modes = allModes{vocabLevelItr-1};
        vocabLevelLabels = [vocabLevel.label];
        
       %% Here, we find edges for the nodes in the vocabulary. 
       % Check the level we need to swith to centroid type edges,
       % no matter what the previous option is.
        if edgeChangeLevel == (vocabLevelItr - 1)
           options.edgeType = 'centroid';
        end
        [allEdges, allEdgeProbs] = extractEdgesInference(nodes, modes, leafNodeArr, firstLevelAdjNodes, options, vocabLevelItr-1);
        if vocabLevelItr == 2
           nonemptyAdjNodeIdx = cellfun(@(x) ~isempty(x), allEdges);
           firstLevelAdjNodes = cell(size(allEdges));
           firstLevelAdjNodes(nonemptyAdjNodeIdx) = cellfun(@(x) x(:,1), allEdges(nonemptyAdjNodeIdx), 'UniformOutput', false);
        end
        
        %% Starting inference from layer 2. 
         for vocabItr = 1:numel(vocabLevel)
             adaptiveThreshold = single(optimalThresholds(vocabLevelItr) * ((size(vocabLevel(vocabItr).adjInfo,1)+1)*2-1)) + singlePrecision; % Hard threshold for cost of matching two subs.
             %% Subgraph matching.
             % Start with the center.
             centerId = vocabLevel(vocabItr).children(1);
             centerMatchCosts = nodeDistanceMatrix(centerId, nodes(:,1))';
             validInstances = centerMatchCosts < adaptiveThreshold;
             
             if nnz(validInstances) == 0
                 continue;
             end

            % Allocate space to hold the instances.
             instanceChildren = int32(find(validInstances));
             instanceMatchCosts = centerMatchCosts(validInstances);
             instancePosProbs = ones(size(instanceChildren,1), 1, 'single');

            %% Get descriptors for edges in the vocabulary node.
             vocabEdges = vocabLevel(vocabItr).adjInfo;

            %% Iteratively, try to match each edge to the instances.
             for edgeItr = 1:size(vocabEdges,1)
                if isempty(instanceChildren) 
                    break;
                end
                childInstances = cell(size(instanceChildren,1),1);
                childMatchingCosts = cell(size(instanceChildren,1),1);
                childPosProbs = cell(size(instanceChildren,1),1);
                for parentItr = 1:size(instanceChildren,1)
                    % Retrieve neighbors close to the center node in
                    % question. We're trying to match an edge to the part
                    % being considered.
                    edges = allEdges{instanceChildren(parentItr)};
                    edgeProbs = allEdgeProbs{instanceChildren(parentItr)};
                    if isempty(edges)
                       continue; 
                    end
                    
                    % Calculate match costs.
                    matchCosts = instanceMatchCosts(parentItr) + edgeDistanceMatrix(edges(:,2), vocabEdges(edgeItr,3)) + ...
                                                                        nodeDistanceMatrix(nodes(edges(:,1),1), ...
                                                                        vocabLevel(vocabItr).children(vocabEdges(edgeItr,2)));                                                                    
                    % If any children are found, we save them.
                    validChildrenIdx = matchCosts < adaptiveThreshold;
                    selectedEdgeProbs = edgeProbs(validChildrenIdx);
                    validChildren = edges(validChildrenIdx,1);
                    
                    if ~isempty(validChildren)
                        newInstanceChildren = zeros(numel(validChildren), size(instanceChildren,2)+1, 'int32');
                        newInstanceChildren(:,1:size(instanceChildren,2)) = repmat(instanceChildren(parentItr,:), numel(validChildren),1);
                        newInstanceChildren(:,end) = validChildren;
                        childInstances(parentItr) = {newInstanceChildren};
                        childMatchingCosts(parentItr) = {matchCosts(validChildrenIdx)};
                        
                        % Assign position probabilities.
                        posProbs = zeros(numel(validChildren), size(instancePosProbs,2) + 1, 'single');
                        posProbs(:,1:(end-1)) = repmat(instancePosProbs(parentItr,:), numel(validChildren),1);
                        posProbs(:,end) = selectedEdgeProbs;
                    else
                         posProbs = [];
                    end
                    childPosProbs(parentItr) = {posProbs};
                end
                instanceChildren = cat(1, childInstances{:});
                instanceMatchCosts = cat(1, childMatchingCosts{:});
                instancePosProbs = cat(1, childPosProbs{:});
                
                % Delete instances with repetitive children.
                numberOfChildren = size(instanceChildren,1);
                validChildren = ones(numberOfChildren,1) > 0;
                for instanceItr = 1:size(instanceChildren,1)
                    if numel(instanceChildren(instanceItr,:)) > ...
                            numel(unique(instanceChildren(instanceItr,:))) 
                        validChildren(instanceItr) = 0;
                    end
                end
                instanceChildren = instanceChildren(validChildren,:);
                instanceMatchCosts = instanceMatchCosts(validChildren,:);
                instancePosProbs = instancePosProbs(validChildren,:);
             end
           
          %% In case of single node subs, we eliminate instances which have outgoing edges.
           if isempty(vocabEdges) && ~isempty(instanceChildren)
               validInstances = cellfun(@(x) isempty(x), allEdges(instanceChildren));
               
               % Eliminate instances which have edges leading out.
               if nnz(validInstances) > 0    
                   instanceChildren = instanceChildren(validInstances, :);
               else
                   instanceChildren = [];
               end
           end

           % Save the instances.
           vocabRealizations(vocabItr) = {instanceChildren};
           
           % Calculating activations here.
           if ~isempty(instanceChildren)
               instanceActivations = ones(size(instanceChildren,1), 1, 'single');
               vocabRealizationsActivations(vocabItr) = {instanceActivations};
           end
        end

        % Get activations.
        activationArr = cat(1, vocabRealizationsActivations{:});
        
        % TODO: Remove this line.
%        activationArr = ones(size(activationArr), 'single');
        numberOfInstances =numel(activationArr);
        
        if numberOfInstances == 0
            break;
        end
        
        %% Create vocabulary realizations.
        maxSize = max(cellfun(@(x) size(x,2), vocabRealizations));
        vocabRealizationsChildren = cellfun(@(x) ...
             cat(2, x, zeros(size(x,1), maxSize - size(x,2), 'int32')), ...
             vocabRealizations, 'UniformOutput', false);
         vocabRealizationsChildren = cat(1, vocabRealizationsChildren{:});
        
        %% We'll form a new nodes array. 
        % Estimate a new activation for each realization by propogating
        % from previous layers.
        newNodes = zeros(numberOfInstances,3);
        newPrecisePositions = zeros(numberOfInstances,2, 'single');
        startIdx = 1;
        for vocabNodeItr = 1:numel(vocabLevel)
           children  = vocabRealizations{vocabNodeItr};
           if isempty(children)
               continue;
           end
           for instanceItr = 1:size(children,1)
               newNodes(startIdx,1) = vocabNodeItr;
               newPrecisePositions(startIdx,:) = precisePositions(children(instanceItr,1),:);
               newNodes(startIdx,2:3) = nodes(children(instanceItr,1), 2:3);
               startIdx = startIdx+1;
           end
        end
        
         %% Assign OR node labels to parts.
         newNodes(:,1) = vocabLevelLabels(newNodes(:,1));
        
         %% Order the nodes in newNodes with real labels first.
         arrayToSort = [newNodes, double(vocabRealizationsChildren)];
         [~, sortedIdx] = sortrows(arrayToSort);
         newNodes = newNodes(sortedIdx, :);
         newPrecisePositions = newPrecisePositions(sortedIdx,:);
         activationArr = activationArr(sortedIdx,:);
         vocabRealizationsChildren = vocabRealizationsChildren(sortedIdx,:);
        
         %% Perform pooling.
         combinedArr = double(newNodes);
         % Sort combinedArr so that it is sorted by decreasing activations.
         [~, idx] = sort(activationArr, 'descend');
         combinedArr = combinedArr(idx,:);
   
         % Downsample the coordinates (pooling), and then perform max operation.
         combinedArr(:,2:3) = floor((combinedArr(:,2:3) - 1)/poolDim) + 1;
         [~, IA, ~] = unique(combinedArr, 'rows', 'stable');

         % Save real indices and activations.
         idx = idx(IA);
         idx = sort(idx);
         activationArr = activationArr(idx);
         vocabRealizationsChildren =vocabRealizationsChildren(idx,:);
         newPrecisePositions = newPrecisePositions(idx,:);
         
         % Create final graphLevel.
         newNodes = newNodes(idx, :);
         newNodes(:,2:3) = int32(floor((double(newNodes(:,2:3)) - 1)/poolDim) + 1);
       
        % Create leaf node array.
        numberOfInstances = size(newNodes,1);
        newLeafNodeArr = cell(numberOfInstances,1);
        for instanceItr = 1:numberOfInstances
             children = vocabRealizationsChildren(instanceItr,:);
             children = children(children>0);
             newLeafNodeArr(instanceItr) = {unique(cat(2, leafNodeArr{children}))};
        end

        % Order nodes.
        [~, sortIdx] = sortrows(newNodes);

        %% Write to output and move on to the next level.
        allNodes(vocabLevelItr) = {newNodes(sortIdx,:)};
        leafNodeArr = newLeafNodeArr(sortIdx);
        allActivations(vocabLevelItr) = {activationArr(sortIdx)};
        nodes = newNodes(sortIdx,:);
        precisePositions = newPrecisePositions(sortIdx,:);
        allPrecisePositions{vocabLevelItr} = precisePositions; 
    end
    activationArr = cat(1, allActivations{:});
    numberOfInstances = sum(cellfun(@(x) size(x,1), allNodes));
    allPrecisePositions = cat(1, allPrecisePositions{:});
    %% If no instances have been found, exit.
    if numberOfInstances<1
        exportArr = [];
        return;
    end
    
    %% Save all nodes to an array, and exit.
    exportArr = zeros(numberOfInstances, 5, 'int32');
    exportArr(:,5) = 1;
    startIdx = 1;
    for vocabLevelItr = 1:numel(allNodes)
        thisLevelNodes = allNodes{vocabLevelItr};
        if isempty(thisLevelNodes)
            break;
        end
        exportArr(startIdx:(startIdx + size(thisLevelNodes,1) - 1), 4) = vocabLevelItr;
        exportArr(startIdx:(startIdx + size(thisLevelNodes,1) - 1), 1:3) = thisLevelNodes(:,1:3);
        startIdx = startIdx + size(thisLevelNodes,1);
    end
end

