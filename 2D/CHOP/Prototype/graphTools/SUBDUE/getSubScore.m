%> Name: getSubScore
%>
%> Description: This function estimates the value of a substructure based on 
%> the metric defined by $evalMetric$. If evalMetric is 'freq', subScore is 
%> simply assigned # of positive instances - # of negative instances. If
%> it is 'size', subScore = (# of positive instances - # of negative
%> instances) * size(sub). If evalMetric is 'mdl', we estimate the compression 
%> amount of an object graph given the list of instance nodes and edges to 
%> compress. It is assumed that the part of the object graphs which are not 
%> directly linked to any of the instances do not change. This criterion makes 
%> it easy to estimate a compression score for the given list of instances, 
%> all of which is considered to belong to a single substructure (part). 
%>
%> @param sub Substructure which will be evaluated.
%> @param allEdges List of all edges in the graph.
%> @param allEdgeNodePairs List of all edge nodes in the graph. 
%> @param allSigns List of all signs of the nodes in the graph.
%> @param mdlNodeWeight Node weight in DL calculations (MDL).
%> @param mdlEdgeWeight Edge weight in DL calculations (MDL).
%> @param overlap If true, overlapping instances are all taken into
%> account. Otherwise only unique instances (in terms of node sets) are
%> considered.
%> @param isMDLExact If true, exact MDL calculation. Approximate otherwise.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.09.2014
function [subScore, sub] = getSubScore(sub, allEdgeCounts, allEdgeNodePairs, evalMetric, ...
                allSigns, mdlNodeWeight, mdlEdgeWeight, overlap, isMDLExact, avgDegree, singleNodeSubThreshold)
    isMultiNodeSub = ~isempty(sub.edges);
    % Read signs and edges of the instances.
    instanceSigns = allSigns(sub.instanceCenterIdx);
    
    % If we're doing something other than MDL: 'size' (sub. size x frequency)
    % or frequency based evaluation, we do not need to proceed with the
    % following steps under following circumstances:
    if ~isMultiNodeSub || ... % Single node sub
            (isMultiNodeSub && ~overlap) || ... % multi node + overlap
            (isMultiNodeSub && strcmp(evalMetric, 'mdl')) % multi node + mdl
        centerIdx = sub.instanceCenterIdx;
        instanceEdgeCounts = allEdgeCounts(centerIdx);
        
        % Calculate outgoing nodes (destinations of edges where
        % instance's children are the source).
        instanceChildren = sub.instanceChildren;
        numberOfNodes = numel(allEdgeCounts);
        numberOfInstances = numel(instanceSigns);   % Beware: Changes if overlap not allowed!

        % If this is a single sub, we only take instances with no POSSIBLE
        % edges into account. If we do not do that, they end up clogging
        % top spots in bestSubs as they usually cover a lot more area than
        % their children.
        if ~isMultiNodeSub
            % Check if this sub qualifies as a single node sub.
            tempIdx = instanceEdgeCounts == 0;
            if nnz(tempIdx) >= numberOfInstances * singleNodeSubThreshold
                 validInstances = ones(numberOfInstances,1) > 0;
            else
                 validInstances = zeros(numberOfInstances,1) > 0;
            end
            instanceSigns = instanceSigns(validInstances);
            instanceChildren = instanceChildren(validInstances);
            
            % Assign valid instances.
            sub.instanceCenterIdx = sub.instanceCenterIdx(validInstances,:);
            sub.instanceChildren = sub.instanceChildren(validInstances,:);
            sub.instanceSigns = sub.instanceSigns(validInstances,:);
        end

        % In multi-node case, if overlaps are not allowed, filter out overlapping instances.
        if isMultiNodeSub && ~overlap
            % Find overlapping instances and remove them by marking as invalid.
            matchedNodes = zeros(numberOfNodes,1)>0;
            validInstances = zeros(numberOfInstances,1) > 0;
            for instItr = 1:numberOfInstances
                children = instanceChildren(instItr,:);
                if nnz(matchedNodes(instanceChildren(instItr,:))) == 0
                    matchedNodes(children) = 1;
                    validInstances(instItr) = 1;
                end
            end

            % Filter out data for invalid instances from existing data structures.
            instanceSigns = instanceSigns(validInstances);
            instanceChildren = instanceChildren(validInstances, :);
        end
    end
    
    % Set constants for each instance. 1 if positive, -1 if negative.
    instanceConstants = ones(numel(instanceSigns),1);
    instanceConstants(~instanceSigns) = -1;
    
    % Calculate value of the substructure.
    switch evalMetric
        case 'freq'
            subScore = sum(instanceConstants);
        case 'size'
            subScore = sum(instanceConstants) * (size(sub.edges,1) + 1);
        case 'mdl'
            if isempty(instanceChildren)
               subScore = 0;
            else
              % Get average degree of children.
                if overlap
                    uniqueNodes = fastsortedunique(sort(instanceChildren));
                else
                    uniqueNodes = sort(instanceChildren(:), 'ascend');
                end

                %% DL calculation for compressed object graph takes place.
                % We only take the modified parts of the new graph into account. It
                % makes things easier.
                % 0) Calculate multiplication constants on each node based on sign. (1
                % if positive, -1 if negative)
                multConstants = ones(numel(allSigns),1);
                multConstants(~allSigns) = -1;

                % 1) Remove unique instance nodes.
                subScore = sum(multConstants(uniqueNodes)) * mdlNodeWeight;

                % 2) Remove unique edges connecting instance nodes to external nodes
                % and in between instances of this sub.
                deletedEdges = ismembc(allEdgeNodePairs(:,1), uniqueNodes) | ismember(allEdgeNodePairs(:,2), uniqueNodes);
                subScore = subScore + sum(multConstants(allEdgeNodePairs(deletedEdges,1))) * mdlEdgeWeight;

                % 3) Add new nodes, one per each new instance.
                subScore = subScore - sum(instanceConstants) * mdlNodeWeight;

                % 4) Add new edges connecting new nodes. In case of exact MDL (actually it is not
                % exact in the sense that the new graph will be passed to the new level
                % as is), a better estimation for the value is used. If not exact, it
                % is assumed that the degree of the graph will be restricted and
                % preserved, so new number of edges will match existing characteristics
                % of the graph. An approximate number of new edges will be added in
                % that case.
                if isMDLExact
                    %% More thorough MDL calculation
                    % Here, put external edges (going in and out new instance
                    % vertices from/to external vertices). No overlap edges are
                    % considered in this DL description. We only link external nodes to
                    % new instance nodes, and no links between new nodes are
                    % considered. The reason to omit this is its immense complexity. We
                    % leave this point as future work. We're open to contributions
                    % towards a better, more correct compression algorithm :)
                    if ~overlap
                        % If overlaps are not allowed, and there are no
                        % background instances, the exact mdl calculation can
                        % be greatly simplified, resulting in a faster process.
                        instChildrenIds = 1:numberOfNodes;
                        for instItr = 1:size(instanceChildren, 1)
                            instChildrenIds(instanceChildren(instItr, :)) = instItr;
                        end
                        instanceNeighbors = allEdgeNodePairs(xor(ismembc(allEdgeNodePairs(:,1), uniqueNodes), ismember(allEdgeNodePairs(:,2),uniqueNodes)),:);
                        instanceNeighbors = instChildrenIds(instanceNeighbors);
                        instanceNeighbors = unique(instanceNeighbors, 'rows');
                        subScore = subScore - round(sum(multConstants(instanceNeighbors(:,1))) * mdlEdgeWeight);
                    else
                        error('Exact MDL case with overlap has not been handled.');
                        % TODO: Will be fixed, currently not working.
                        instanceOutNeighbors = cellfun(@(x) numel(setdiff(allEdgeNodePairs(ismember(allEdgeNodePairs(:,1), x),2), x)), instanceChildren);
                        instanceInNeighbors = cellfun(@(x) numel(setdiff(allEdgeNodePairs(ismember(allEdgeNodePairs(:,2), x),1), x)), instanceChildren);
                        subScore = subScore - round(sum((instanceOutNeighbors + instanceInNeighbors) .* instanceConstants) * mdlEdgeWeight);
                    end
                else
                    %% Approximate MDL calculation
                    % Adding a number of edges per each new instance, which is
                    % consistent with the average degree of the graph. The basic
                    % assumption is that the degree/topology of the graph will not
                    % change. The number of in/out edges from/to external edges is
                    % assumed to be the same, that's why the average degree is
                    % multiplied by two.
                    subScore = subScore - round(2 * avgDegree * sum(instanceConstants) * mdlEdgeWeight);
                end

                % 5) Finally, consider the own description length of the sub.
                subScore = subScore - ((size(sub.edges,1) + 1) * mdlNodeWeight + ...
                                                       size(sub.edges,1) * mdlEdgeWeight);
            end
        otherwise
            display('Error: evalMetric not implemented (in getSubScore.m).');
    end
end