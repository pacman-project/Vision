%> Name: extendSub
%>
%> Description: extendSub(..) extends 'sub' in all possible ways
%> by extending its instances, and returns the new sub list 'extendedSubs'
%> along with instances of each sub in returned list.
%> 
%> @param sub Sub to be extended.
%> @param allEdges List of all edges in the graph (Nx4).
%>
%> @retval extendedSubs Extended sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function [extendedSubs] = extendSub(sub, allEdges, nodeDistanceMatrix, edgeDistanceMatrix, threshold)
    
    centerIdx = sub.instanceCenterIdx;
    subAllEdges = {allEdges(centerIdx).adjInfo}';
    % Get unused edges from sub's instances.
    allUsedEdgeIdx = sub.instanceEdges;
    if ~isempty(allUsedEdgeIdx)
        allUsedEdgeIdx = mat2cell(allUsedEdgeIdx, ones(size(allUsedEdgeIdx,1),1), size(allUsedEdgeIdx,2));
        allUnusedEdgeIdx = cellfun(@(x,y) fast_setdiff(1:size(x,1),y), subAllEdges, allUsedEdgeIdx, 'UniformOutput', false);
        allUnusedEdges = cellfun(@(x,y) x(y,:), subAllEdges, allUnusedEdgeIdx, 'UniformOutput', false);
        allUnusedEdgeIdx = [allUnusedEdgeIdx{:}]';
    else
        allUnusedEdges = subAllEdges;
        allUnusedEdgeIdx = [];
    end
    
    % Record which edge belongs to which instance. 
    allEdgeInstanceIds = zeros(sum(cellfun(@(x) size(x,1), allUnusedEdges)),1);
    allEdgePrevCosts = zeros(size(allEdgeInstanceIds), 'single');
    allEdgeExactMatchFlags = zeros(size(allEdgeInstanceIds), 'uint8')>0;
    itrOffset = 1;
    unusedEdgeCounts = cellfun(@(x) size(x,1), allUnusedEdges);
    for itr = 1:numel(allUnusedEdges)
        beginOffset = itrOffset;
        endOffset = (beginOffset+(unusedEdgeCounts(itr)-1));
        allEdgeInstanceIds(beginOffset:endOffset) = itr;
        allEdgePrevCosts(beginOffset:endOffset) = sub.instanceMatchCosts(itr);
        allEdgeExactMatchFlags(beginOffset:endOffset) = sub.instanceExactMatchFlags(itr);
        itrOffset = itrOffset + unusedEdgeCounts(itr);
    end         
    allUnusedEdges = cat(1, allUnusedEdges{:});
    
    % If no edges remain, exit.
    if isempty(allUnusedEdges)
        extendedSubs = [];
        return; 
    end
    
    % Eliminate the edges which exist only in validation data. We do not
    % enumerate any edges which do not exist in training data.
    enumeratedEdges = allUnusedEdges(allEdgeExactMatchFlags, 3:4);
    
    % Get unique rows of [edgeLabel, secondVertexLabel]
    uniqueEdgeTypes = unique(enumeratedEdges, 'rows');
    
    % Eliminate roughly half of the unique edge types but lexicographical
    % ordering. If new edge is lexicographically lower, it is not
    % enumerated.
    edges = sub.edges;
    if ~isempty(edges)
       uniqueEdgeTypesPart1 = uniqueEdgeTypes(uniqueEdgeTypes(:,1) == max(edges(:,1)), :);
       if ~isempty(uniqueEdgeTypesPart1)
           maxNode = max(edges(edges(:,1) == max(edges(:,1)),2));
           uniqueEdgeTypesPart1 = uniqueEdgeTypesPart1(uniqueEdgeTypesPart1(:,2) > maxNode, :);
       end
       uniqueEdgeTypesPart2 = uniqueEdgeTypes(uniqueEdgeTypes(:,1) > max(edges(:,1)), :);
       uniqueEdgeTypes = cat(1, uniqueEdgeTypesPart1, uniqueEdgeTypesPart2);
    end
    
    % Extend the definition in sub with each edge type in uniqueEdgeTypes.
    % In addition, we pick suitable instances, add this edge, and mark used
    % field of relevant instances.
    numberOfEdgeTypes = size(uniqueEdgeTypes,1);
    if numberOfEdgeTypes == 0
        extendedSubs = [];
        return;
    end
    
    % Save local indices for all edges, to be used later.
    allLocalEdgeIdx = cellfun(@(x) 1:size(x,1), subAllEdges, 'UniformOutput', false);
    allLocalEdgeIdx = [allLocalEdgeIdx{:}]';
    
    % Allocate space for new subs and fill them in.
    validSubs = ones(numberOfEdgeTypes,1) > 0;
    extendedSubs(numberOfEdgeTypes) = Substructure;
    for edgeTypeItr = 1:numberOfEdgeTypes
        % Assign sub-definition type and other info
        newSub = Substructure;
        newSub.edges = sub.edges;
        newSub.centerId = sub.centerId;
        newSub.mdlScore = 0;
        newSub.edges = [newSub.edges; uniqueEdgeTypes(edgeTypeItr,:)];
        
        %% Find instances of this new sub, and mark the new edges as 'used' in each of its subs.
        % Extend the subs by taking the threshold into account. Unless the
        % instance's collective matching cost surpasses the threshold, it
        % is a valid instance.
        newEdgeCosts = edgeDistanceMatrix(uniqueEdgeTypes(edgeTypeItr,1), allUnusedEdges(:,3))';
        newNodeCosts = nodeDistanceMatrix(uniqueEdgeTypes(edgeTypeItr,2), allUnusedEdges(:,4))';
        edgesToExtendCosts = allEdgePrevCosts + newEdgeCosts + newNodeCosts;
        edgesToExtendIdx = edgesToExtendCosts < threshold;
        
        % Check for exact matching again.
        exactMatchFlags = allEdgeExactMatchFlags & ...
             min(edgeDistanceMatrix(uniqueEdgeTypes(edgeTypeItr,1),:)) == newEdgeCosts & ...
             min(nodeDistanceMatrix(uniqueEdgeTypes(edgeTypeItr,2),:)) == newNodeCosts;
        
        % Save instance ids.
        edgeInstanceIds = allEdgeInstanceIds(edgesToExtendIdx);
        allChildren = [sub.instanceChildren(edgeInstanceIds,:), ...
             allUnusedEdges(edgesToExtendIdx,2)];
%        [allChildren, allChildrenSortIdx] = sort(allChildren, 2);

        % Calculate mappings of instances' nodes to the description.
        newMappings = [sub.instanceMappings(edgeInstanceIds, :), ...
            ones(size(edgeInstanceIds,1), 1, 'uint8') * size(allChildren,2)];
%         [rows, cols] = size(allChildren);
%         R = repmat((1:rows)',[1 cols]);
%         nIdx = R + (allChildrenSortIdx-1)*rows;
%         newMappings = newMappings(nIdx);
        
        %% Note: allChildren may have duplicate entries, since an instance can
        % have multiple parse trees. We handle these cases by only keeping
        % unique instances. In addition, for each instance, the minimum
        % cost of matching is kept here.
        curInstanceMatchCosts = edgesToExtendCosts(edgesToExtendIdx);
        curInstanceExactMatchFlags = exactMatchFlags(edgesToExtendIdx);
        [minMatchCosts, sortIdx] = sort(curInstanceMatchCosts, 'ascend');
        curInstanceExactMatchFlags = curInstanceExactMatchFlags(sortIdx);
        sortedAllChildren = allChildren(sortIdx, :);
        sortedCenterIdx = sub.instanceCenterIdx(edgeInstanceIds(sortIdx));

        % Eliminate instances which have the same node set, and the same
        % center node. 
 %       [~, validIdx, ~] = unique([sortedCenterIdx, sortedAllChildren], 'rows', 'stable');
        validIdx = (1:numel(sortedCenterIdx))';
        
        % Keep ordering but remove duplicates.
        sortIdx = sortIdx(validIdx);
        
        % Get minimum matching costs and children.
        minMatchCosts = minMatchCosts(validIdx, :);
        curInstanceExactMatchFlags = curInstanceExactMatchFlags(validIdx,:);
        sortedAllChildren = sortedAllChildren(validIdx, :);
        
        % Finally, order children by rows.
        [newSub.instanceChildren, idx] = sortrows(sortedAllChildren);
        sortIdx = sortIdx(idx);
        edgeInstanceIds = edgeInstanceIds(sortIdx, :);
        newMappings = newMappings(sortIdx, :);
        
        %% Assign all relevant instance-related fields of the sub.
        newSub.instanceMappings = newMappings;
        newSub.instanceCenterIdx = sub.instanceCenterIdx(edgeInstanceIds);
        newSub.instanceSigns = sub.instanceSigns(edgeInstanceIds);
        newSub.instanceCategories = sub.instanceCategories(edgeInstanceIds);
        newSub.instanceMatchCosts = minMatchCosts(idx,:);
        newSub.instanceExactMatchFlags = curInstanceExactMatchFlags(idx,:);
        newSub.instanceValidationIdx = sub.instanceValidationIdx(edgeInstanceIds);
        
        % Add the edge to the definition.
        existingEdges = sub.instanceEdges;
        combinedEdges = zeros(numel(edgeInstanceIds), size(existingEdges,2)+1, 'uint8');
        if ~isempty(existingEdges)
            existingEdges = existingEdges(edgeInstanceIds,:);
            combinedEdges(:, 1:size(existingEdges,2)) = existingEdges;
        end
        
        % Adding the local edges to the instance edge lists.
        edgesToExtendLinIdx = find(edgesToExtendIdx);
        if ~isempty(sub.edges)
            relevantLocalEdgeIdx = allUnusedEdgeIdx(edgesToExtendLinIdx(sortIdx));
        else
            relevantLocalEdgeIdx = allLocalEdgeIdx(edgesToExtendLinIdx(sortIdx));
        end
        combinedEdges(:,end) = relevantLocalEdgeIdx;
        newSub.instanceEdges = combinedEdges;
        
        % All instances assigned, good to go.
        extendedSubs(edgeTypeItr) = newSub;
    end
    
    extendedSubs = extendedSubs(validSubs);
end

function Z = fast_setdiff(X,Y)
       check = false(1, max(max(X), max(Y)));
       check(X) = true;
       check(Y) = false;
       Z = X(check(X));  
end