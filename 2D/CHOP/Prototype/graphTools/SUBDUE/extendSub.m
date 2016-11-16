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
function [extendedSubs] = extendSub(sub, allEdges, allEdgeCounts, singleInstanceFlag)
    % Get center list.
    centerIdx = sub.instanceCenterIdx;
    
    % Record which edge belongs to which instance. 
    edgeCounts = allEdgeCounts(centerIdx);
    allEdgeInstanceIds = zeros(sum(edgeCounts),1);
    itrOffset = 1;
    for itr = 1:numel(centerIdx)
        beginOffset = itrOffset;
        allEdgeInstanceIds(beginOffset:(beginOffset+(edgeCounts(itr)-1))) = itr;
        itrOffset = itrOffset + edgeCounts(itr);
    end         
    subAllEdges = cat(1, allEdges{centerIdx});
    
    % If there are no edges, exit.
    if isempty(subAllEdges)
        extendedSubs = [];
        return; 
    end
    
    % Eliminate the edges which exist only in validation data. We do not
    % enumerate any edges which do not exist in training data.
    enumeratedEdges = subAllEdges(:, 3:4);
    
    % Get unique rows of [edgeLabel, secondVertexLabel]
    uniqueEdgeTypes = unique(enumeratedEdges, 'rows');
    
    instanceChildren = sub.instanceChildren;
    
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
 %      uniqueEdgeTypesPart2 = uniqueEdgeTypes(uniqueEdgeTypes(:,1) >= max(edges(:,1)), :);
       uniqueEdgeTypes = cat(1, uniqueEdgeTypesPart1, uniqueEdgeTypesPart2);
    end
    
    % If no edges remain, exit.
    if isempty(uniqueEdgeTypes)
       extendedSubs = [];
       return;
    end
    
    % Extend the definition in sub with each edge type in uniqueEdgeTypes.
    % In addition, we pick suitable instances, add this edge, and mark used
    % field of relevant instances.
    numberOfEdgeTypes = size(uniqueEdgeTypes,1);
    
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
        edgesToExtendIdx = subAllEdges(:,3) == uniqueEdgeTypes(edgeTypeItr,1) & ...
            subAllEdges(:,4) == uniqueEdgeTypes(edgeTypeItr,2);
        
        % If single instance subs are not allowed, we render this sub
        % invalid if it has a single instance.
        if ~singleInstanceFlag && nnz(edgesToExtendIdx) == 1
             validSubs(edgeTypeItr) = 0;
             continue;
        end
        
        % Save instance ids.
        edgeInstanceIds = allEdgeInstanceIds(edgesToExtendIdx);
        allChildren = [instanceChildren(edgeInstanceIds,:), ...
             subAllEdges(edgesToExtendIdx,2)];

        % Finally, order children by rows.
        [sortedChildren, idx] = sortrows(allChildren);
        
%         % Go through the children and select only cliques (all nodes
%         % connected)
%         edgeCount = size(sortedChildren,2) - 1;
%         if edgeCount > 1
%               seedRow = sortedChildren(:,end);
%               for edgeItr = 1:(edgeCount-1)
%                    comparedRow = sortedChildren(:,1+edgeItr);
%                    tempIdx = seedRow + (comparedRow-1) * nodeCount;
%                    if edgeItr == 1
%                         validIdx = full(adjMatrix(tempIdx));
%                    else
%                         validIdx = validIdx & full(adjMatrix(tempIdx));
%                    end
%               end
%               idx = idx(validIdx);
%         end
%         
%         if isempty(idx)
%              validSubs(edgeTypeItr) = 0;
%              continue;
%         end
        
        % Keep valid instances.
        newSub.instanceChildren= sortedChildren;
        edgeInstanceIds = edgeInstanceIds(idx, :);
        
        %% Assign all relevant instance-related fields of the sub.
        newSub.instanceCenterIdx = sub.instanceCenterIdx(edgeInstanceIds);
        newSub.instanceSigns = sub.instanceSigns(edgeInstanceIds);
        
        % All instances assigned, good to go.
        extendedSubs(edgeTypeItr) = newSub;
    end
    
    extendedSubs = extendedSubs(validSubs);
end