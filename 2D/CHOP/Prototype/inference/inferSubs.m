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
function [graphLevel] = inferSubs(vocabLevel, graphLevel, nodeDistanceMatrix, edgeDistanceMatrix, threshold)
    % Read data into helper data structures.
    edges = cat(1, graphLevel.adjInfo);
    
    % If no edges are present, time to return.
    if isempty(edges)
        graphLevel = [];
        return;
    end
    
    labelIds = cat(1, graphLevel.labelId);
    
    %% Match subs from vocabLevel to their instance in graphLevel.
    vocabRealizations = cell(numel(vocabLevel),1);
    for vocabItr = 1:numel(vocabLevel)
         adaptiveThreshold = single(threshold * ((size(vocabLevel(vocabItr).adjInfo,1)+1)*2-1)) + 0.0001; % Hard threshold for cost of matching two subs.
         %% Subgraph matching.
         % Start with the center.
         centerId = vocabLevel(vocabItr).children(1);
         centerMatchCosts = nodeDistanceMatrix(centerId, labelIds)';
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
                currentEdges = graphLevel(instanceChildren(parentItr,1)).adjInfo;
                if size(instanceChildren,2) > 1
                    currentEdges = currentEdges(~ismember(currentEdges(:,2), instanceChildren(parentItr,2:end)), :);
                end
                currentEdges = [currentEdges(:,3), currentEdges(:,2)];
                matchCosts = instanceMatchCosts(parentItr) + edgeDistanceMatrix(currentEdges(:,1), vocabEdges(edgeItr,3)) + ...
                                                                    nodeDistanceMatrix(labelIds(currentEdges(:,2)), vocabLevel(vocabItr).children(vocabEdges(edgeItr,2)));
                validChildrenIdx = matchCosts < adaptiveThreshold;
                validChildren = currentEdges(validChildrenIdx,2);
                
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
        end
       vocabRealizations(vocabItr) = {instanceChildren};
    end
    numberOfInstanceArr = cellfun(@(x) size(x,1), vocabRealizations);
    numberOfInstances = sum(numberOfInstanceArr);
    clear graphLevel;
    
    %% If no instances have been found, exit.
    if numberOfInstances<1
        graphLevel = [];
        return;
    end
    
    %% Generate graph level.
    graphLevel(numberOfInstances) = GraphNode;
    
    % Collect label ids, children and signs to assign to instances.
    labelIds = zeros(numberOfInstances,1, 'int32');
    instanceOffset = 1;
    for vocabItr = 1:numel(vocabLevel)
       if numberOfInstanceArr(vocabItr) > 0
            instanceEndOffset = instanceOffset + (numberOfInstanceArr(vocabItr)-1);
            labelIds(instanceOffset:instanceEndOffset) = int32(vocabItr);
            instanceOffset = instanceEndOffset+1;
       end
    end
    [~, sortIdx] = sort(labelIds);
    labelIds = num2cell(labelIds);
    allChildren = cellfun(@(x) mat2cell(x, ones(size(x,1),1), size(x,2)), vocabRealizations, 'UniformOutput', false);
    allChildren = cat(1, allChildren{:});
    
    %% Assign labelId, children and sign of the children.
    [graphLevel.labelId] = deal(labelIds{:});
    [graphLevel.children] = deal(allChildren{:});
    [graphLevel.sign] = deal(1);
    
    % Sort graph level based on label ids.
    graphLevel = graphLevel(sortIdx);
    
    % Get unique instances (to be consistent with training)
    allChildren = allChildren(sortIdx);
    allChildren = cellfun(@(x) mat2str(sort(x)), allChildren, 'UniformOutput', false);
    [~, IC, ~] = unique(allChildren, 'stable');
    graphLevel = graphLevel(IC);
end

