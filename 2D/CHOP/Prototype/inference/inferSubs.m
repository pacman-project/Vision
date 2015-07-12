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
function [exportArr, activationArr] = inferSubs(vocabulary, nodes, nodeActivations, distanceMatrices, optimalThresholds, vocabUpdatedLabels, options)
    % Read data into helper data structures.
    edgeDistanceMatrix = double(options.edgeDistanceMatrix);
    noveltyThr = 1 - options.noveltyThr;
    firstLevelAdjNodes = [];
    
    % If fast inference is not required, we do not perform inhibition.
    if ~(options.fastInference)
       noveltyThr = 1; 
    end
    multiplier = round(1/options.singlePrecision);
    singlePrecision = options.singlePrecision;
    
    exportArr = [];
    if isempty(nodes) || isempty(vocabulary)
        return;
    end
    allNodes = cell(numel(vocabulary),1);
    allActivations = cell(numel(vocabulary),1);
    allNodes(1) = {nodes};
    allActivations(1) = {nodeActivations};
    leafNodeArr = num2cell((int32(1:size(nodes,1)))');
    
    for vocabLevelItr = 2:numel(vocabulary)
        prevActivations = allActivations{vocabLevelItr-1};
        scale = (1/options.scaling)^(vocabLevelItr-2);
        neighborhood = fix(options.edgeRadius * scale);
        
        %% Match subs from vocabLevel to their instance in graphLevel.
        vocabLevel = vocabulary{vocabLevelItr};
        vocabRealizations = cell(numel(vocabLevel),1);
        vocabRealizationsActivation = cell(numel(vocabLevel),1);
        nodeDistanceMatrix = distanceMatrices{vocabLevelItr-1};
        
       %% Here, we find edges for the nodes in the vocabulary. 
        allEdges = extractEdgesInference(nodes, leafNodeArr, firstLevelAdjNodes, options, vocabLevelItr-1);
        if vocabLevelItr == 2
           firstLevelAdjNodes = cellfun(@(x) x(:,1), allEdges, 'UniformOutput', false);
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

            %% Get descriptors for edges in the vocabulary node.
             vocabEdges = vocabLevel(vocabItr).adjInfo;

            %% Iteratively, try to match each edge to the instances.
             for edgeItr = 1:size(vocabEdges,1)
                if isempty(instanceChildren) 
                    break;
                end
                childInstances = cell(size(instanceChildren,1),1);
                childMatchingCosts = cell(size(instanceChildren,1),1);
                for parentItr = 1:size(instanceChildren,1)
                    % Retrieve neighbors close to the center node in
                    % question. We're trying to match an edge to the part
                    % being considered.
                    edges = allEdges{instanceChildren(parentItr)};
                    if isempty(edges)
                       continue; 
                    end
                    
                    % Calculate match costs.
                    matchCosts = instanceMatchCosts(parentItr) + edgeDistanceMatrix(edges(:,2), vocabEdges(edgeItr,3)) + ...
                                                                        nodeDistanceMatrix(nodes(edges(:,1),1), ...
                                                                        vocabLevel(vocabItr).children(vocabEdges(edgeItr,2)));                                                                    
                    % If any children are found, we save them.
                    validChildrenIdx = matchCosts < adaptiveThreshold;
                    validChildren = edges(validChildrenIdx,1);

                    if ~isempty(validChildren)
                        newInstanceChildren = zeros(numel(validChildren), size(instanceChildren,2)+1, 'int32');
                        newInstanceChildren(:,1:size(instanceChildren,2)) = repmat(instanceChildren(parentItr,:), numel(validChildren),1);
                        newInstanceChildren(:,end) = validChildren;
                        childInstances(parentItr) = {newInstanceChildren};
                        childMatchingCosts(parentItr) = {matchCosts(validChildrenIdx)};
                    end
                end
                instanceChildren = cat(1, childInstances{:});
                instanceMatchCosts = cat(1, childMatchingCosts{:});
                
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
             end
             instanceChildren = sort(instanceChildren, 2);
             
           %% Eliminating duplicate entries in instanceChildren.
            % We handle these cases by only keeping
            % unique instances. In addition, for each instance, the minimum
            % cost of matching is kept here.
            [minMatchCosts, sortIdx] = sort(instanceMatchCosts, 'ascend');
            sortedChildren = instanceChildren(sortIdx, :);
            [~, validIdx, ~] = unique(sortedChildren, 'rows', 'stable');

            % Get minimum matching costs and children.
            instanceMatchCosts = minMatchCosts(validIdx, :);
            sortedChildren = sortedChildren(validIdx, :);

            % Finally, order children by rows.
            [instanceChildren, idx] = sortrows(sortedChildren);
            instanceMatchCosts = instanceMatchCosts(idx);
           
          %% In case of single node subs, we eliminate instances which have outgoing edges.
           if isempty(vocabEdges) && ~isempty(instanceChildren)
               validInstances = cellfun(@(x) isempty(x), allEdges(instanceChildren));
               
               % Eliminate instances which have edges leading out.
               if nnz(validInstances) > 0    
                   instanceChildren = instanceChildren(validInstances, :);
                   instanceMatchCosts = instanceMatchCosts(validInstances, :);
               else
                   instanceChildren = [];
                   instanceMatchCosts = [];
               end
           end

           % Save the instances.
           vocabRealizations(vocabItr) = {instanceChildren};
           
           % Calculating activations here.
           if ~isempty(instanceChildren)
               instanceMatchScores = fix(multiplier * ((adaptiveThreshold - instanceMatchCosts) / adaptiveThreshold)) / multiplier;
               instancePrevActivations = prevActivations(instanceChildren);
               if numel(instanceMatchScores)>1
                   instanceMeanActivations = mean(instancePrevActivations,2);
               else
                   instanceMeanActivations = mean(instancePrevActivations);
               end
               
               instanceActivations = instanceMatchScores .* instanceMeanActivations;
    %           instanceActivations = instanceMatchCosts;
               vocabRealizationsActivation(vocabItr) = {instanceActivations};
           end
        end

        % Get activations.
        activationArr = cat(1, vocabRealizationsActivation{:});
        numberOfInstances =numel(activationArr);
        
        if numberOfInstances == 0
            break;
        end
        
        %% We'll form a new nodes array. 
        % Estimate a new activation for each realization by propogating
        % from previous layers.
        newNodes = zeros(numberOfInstances,3);
        startIdx = 1;
        for vocabNodeItr = 1:numel(vocabLevel)
           children  = vocabRealizations{vocabNodeItr};
           if isempty(children)
               continue;
           end
           for instanceItr = 1:size(children,1)
               newNodes(startIdx,1) = vocabNodeItr;
               newNodes(startIdx,2:3) = round(mean(nodes(children(instanceItr,:),2:3),1));
               startIdx = startIdx+1;
           end
        end
        
        %% We do two things here in order to eliminate number of nodes.]
        % Step zero: Get the list of leaf nodes for every instance.
        instanceOffset = 1;
        newLeafNodes = cell(numberOfInstances,1);
        for vocabNodeItr = 1:numel(vocabLevel)
            if ~isempty(vocabRealizations{vocabNodeItr})
                children = vocabRealizations{vocabNodeItr};
                for instanceItr = 1:size(children,1)
                    leafNodes = unique(cat(2, leafNodeArr{children(instanceItr,:)}));
                    newLeafNodes{instanceOffset} = leafNodes;
                    instanceOffset = instanceOffset + 1;
                end
            end
        end
         maxSize = max(cellfun(@(x) size(x,2), vocabRealizations));
         vocabRealizationsLabelIds = num2cell((1:numel(vocabLevel)))';
         vocabRealizationsDescriptors = cellfun(@(x, y) ...
             cat(2, ones(size(x,1), 1, 'int32') * y, x, zeros(size(x,1), maxSize - size(x,2), 'int32')), ...
             vocabRealizations, vocabRealizationsLabelIds, 'UniformOutput', false);
         vocabRealizationsDescriptors = cat(1, vocabRealizationsDescriptors{:});
        
        %% Inhibition! This is the inhibition step that is identical to applyTestInhibition procedure.
        % If that part changes, this should change as well.
        % First, we sort nodes in descending order by their activation
        % scores.
        [sortedActivationArr, sortIdx] = sort(activationArr, 'descend');
        sortedLeafNodeArr = newLeafNodes(sortIdx);
        sortedNewNodes = newNodes(sortIdx, :);
        vocabRealizationsDescriptors = vocabRealizationsDescriptors(sortIdx, :);
        
        % Next, perform inhibition.
%        if noveltyThr < 1 && vocabLevelItr < numel(vocabulary)
        if noveltyThr < 1
            imagePreservedNodes = ones(numberOfInstances,1)>0;
            imageNodeCoords = sortedNewNodes(:,2:3);
            maxSharedLeafNodes = cellfun(@(x) numel(x) * noveltyThr , sortedLeafNodeArr, 'UniformOutput', false);
            for nodeItr = 1:(numberOfInstances-1)
                %% If nobody has erased this node before, it has a right to be in the final graph.
                if imagePreservedNodes(nodeItr) == 0
                    continue;
                end

                %% Get each neighboring node.
                thisNodeCoords = sortedNewNodes(nodeItr,2:3);
                centerArr = repmat(thisNodeCoords, numberOfInstances, 1);
                distances = sqrt(sum((centerArr - imageNodeCoords).^2, 2));
                adjacentNodes = imagePreservedNodes & distances <= neighborhood; 
                adjacentNodes(1:nodeItr) = 0;
                selfLeafNodes = sortedLeafNodeArr{nodeItr};

                %% Go over each adjacent node, and apply inhibition if their leaf nodes are too common under current novelty threshold.
                imagePreservedNodes(adjacentNodes) = cellfun(@(x,y) sum(ismembc(x, selfLeafNodes)) <= y, ...
                  sortedLeafNodeArr(adjacentNodes), maxSharedLeafNodes(adjacentNodes));
            end
            sortedNewNodes = sortedNewNodes(imagePreservedNodes, :);
            sortedLeafNodeArr = sortedLeafNodeArr(imagePreservedNodes, :);
            sortedActivationArr = sortedActivationArr(imagePreservedNodes, :);
            vocabRealizationsDescriptors = vocabRealizationsDescriptors(imagePreservedNodes,:);
        end
        
        % Update the graph label ids, based on the new ranking.
        % Compositions are ordered in a different way in the vocabulary. 
        % Relevant info can be found in postProcessParts.m.
        updatedOrdering = vocabUpdatedLabels{vocabLevelItr-1};
        vocabRealizationsDescriptors(:,1) = updatedOrdering(vocabRealizationsDescriptors(:,1));
        sortedNewNodes(:,1) = updatedOrdering(sortedNewNodes(:,1));
        
        % Finally, sort everything back.
        [~, idx] = sortrows(vocabRealizationsDescriptors);
        sortedNewNodes = sortedNewNodes(idx, :);
        sortedLeafNodeArr = sortedLeafNodeArr(idx,:);
        sortedActivationArr = sortedActivationArr(idx,:);
       
        % Write to output and move on to the next level.
        allNodes(vocabLevelItr) = {sortedNewNodes};
        allActivations(vocabLevelItr) = {sortedActivationArr};
        leafNodeArr = sortedLeafNodeArr;
        nodes = sortedNewNodes;
    end
    activationArr = cat(1, allActivations{:});
    
    numberOfInstances = sum(cellfun(@(x) size(x,1), allNodes));
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
        exportArr(startIdx:(startIdx + size(thisLevelNodes,1) - 1), 1:3) = thisLevelNodes;
        startIdx = startIdx + size(thisLevelNodes,1);
    end
end

