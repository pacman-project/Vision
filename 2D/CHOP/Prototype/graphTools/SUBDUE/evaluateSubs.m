%> Name: evaluateSubs
%>
%> Description: The evaluation function of SUBDUE. Based on the
%> evaluation metric, the value of each substructure in subs list is
%> calculated, and saved within subs. The MDL calculation takes place here.
%> Key points:
%>  1) DL estimation not rigorous. It is considered that 
%>        each node needs a label (int) and a pointer to its edges (int)
%>        each edge needs a node label (int) for its destination node, an
%>        edge label(int) and a binary isDirected label (bit)
%>     in the final graph description. Node and edge weights in DL
%>     calculation is stored in options.
%> 
%> @param subs Sub list which will be evaluated.
%> @param evalMetric Evaluation metric.
%> @param allEdges List of all edges in the graph.
%> @param allEdgeNodePairs List of all edge node pairs in the graph.
%> @param allSigns List of all signs of the nodes in the graph.
%> @param graphSize Size of the graph.
%> @param overlap If true, overlapping instances are considered in
%> evaluation of a sub.
%> @param mdlNodeWeight Node weight in DL calculations (MDL).
%> @param mdlEdgeWeight Edge weight in DL calculations (MDL).
%> @param isMDLExact If true, exact MDL calculation. Approximate otherwise.
%>
%> @retval subs Evaluated sub list.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 24.02.2014
%> Ver 1.1 on 01.09.2014 Removal of global parameters.
function [subs, validSubs, validExtSubs] = evaluateSubs(subs, evalMetric, allEdgeCounts, allEdgeNodePairs, allSigns, allCoords, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, ...
     allLeafNodes, level1CoordsPooled, halfRFSize, minRFCoverage, maxShareability, maxLeafCounts, avgDegree, singleNodeSubThreshold)
    numberOfSubs = numel(subs);
    validSubs = ones(numberOfSubs,1) > 0;
    validExtSubs = ones(numberOfSubs,1) > 0;
    maxCover = 0.9;
    for subItr = 1:numberOfSubs
        
        % Calculate Substructure score.
        [subScore, sub] = getSubScore(subs(subItr), allEdgeCounts, allEdgeNodePairs, evalMetric, ...
           allSigns, mdlNodeWeight, mdlEdgeWeight, overlap, isMDLExact, avgDegree, singleNodeSubThreshold);
       
       % If this sub has multiple children, we process the instances to
       % make sure those instances which cover a large portion of their
       % receptive field survives.
        if size(subs(subItr).instanceChildren,2) > 1 && (minRFCoverage > 0 || maxShareability < 1)
              allChildren = subs(subItr).instanceChildren;
              % Here, we check if we're actually covering enough of every
              % receptive field out there.
              numberOfInstances = size(allChildren,1);
              numberOfChildren = size(allChildren,2);
              coveredLeafNodes = cell(numberOfInstances,1);
              shareabilityArr = zeros(numberOfInstances,1);
%               for childItr = 1:numberOfInstances
%                    %% A different operation here checks if added child is suitable for the others.
% %                    if numberOfChildren > 2
% %                         leafNodeLists = allLeafNodes(allChildren(childItr,2:end));
% %                         comparedList = leafNodeLists{end};
% %                         for itr = 1:(numberOfChildren-2)
% %                              curList = leafNodeLists{itr};
% %                              if numel(fastintersect(comparedList, curList)) > min(maxShareability * numel(comparedList), maxShareability * numel(curList))
% %                                   shareabilityArr(childItr) = 1;
% %                                   break;
% %                              end
% %                         end
% %                    end
%                    
% %                    %% Calculate support, and determine if we should eliminate this instance.
% %                    % We compare peripheral node supports to decide whether
% %                    % to eliminate this instance or not. 
% %                    if maxShareability < 1 && numberOfChildren > 2
% %                         numberOfUniqueLeafNodes = numel(tempLeafNodes);
% %                         numberOfRepetitions = numel(tempArr) - numel(tempLeafNodes);
% % 
% %                         % Determine if this instance is any good.
% %                         repetitionRatio = numberOfRepetitions / numberOfUniqueLeafNodes;
% %                         if repetitionRatio > maxShareability * (numberOfChildren-1)
% %                              shareabilityArr(childItr) = 1;
% %                         elseif repetitionRatio > maxShareability
% %                              tempArr = [false, diff(tempArr)== 0];
% %                              tempArr = find(tempArr);
% %                              repeatedNodes = numel(tempArr);
% %                              for itr = 2:numel(tempArr)
% %                                    if tempArr(itr) == tempArr(itr-1)+1
% %                                         repeatedNodes = repeatedNodes - 1;
% %                                    end
% %                              end
% %                              shareabilityArr(childItr) = repeatedNodes / numberOfUniqueLeafNodes;
% %                          end
% %                    end
%                    
%                     if minRFCoverage > 0

%                     end
%               end
              
              %% Mark valid instances.
%              if minRFCoverage > 0
%                     tempArr = sort(cat(2, allLeafNodes{allChildren(childItr,:)}));
%                     tempLeafNodes = fastsortedunique(tempArr);
%                     tempLeafNodeCoords = level1CoordsPooled(tempLeafNodes,:);
%                     nodeCoords = allCoords(allChildren(childItr,1),:);
%       %                    tempLeafNodes2 = tempLeafNodes(pdist2(single(tempLeafNodeCoords), single(nodeCoords)) < halfRFSize);
%                     tempLeafNodes = tempLeafNodes(tempLeafNodeCoords(:,1) > nodeCoords(1) - halfRFSize & ...
%                     tempLeafNodeCoords(:,1) < nodeCoords(1) + halfRFSize & ...
%                     tempLeafNodeCoords(:,2) > nodeCoords(2) - halfRFSize & ...
%                     tempLeafNodeCoords(:,2) < nodeCoords(2) + halfRFSize);
%                     coveredLeafNodes{childItr} = int32(tempLeafNodes);
%                   % Find the intersection of two sets, to assess coverage.
%                   coveredLeafNodeCount = cellfun(@(x) numel(x), coveredLeafNodes);
%                   maxCoverLeafNodeCount = maxLeafCounts(allChildren(:,1)); %#ok<PFBNS>
%                   coverageRatios =  coveredLeafNodeCount ./ maxCoverLeafNodeCount;
%                   validInstances = coverageRatios >= minRFCoverage & shareabilityArr <= maxShareability;
%                   
%                    % If the coverage ratios are too high, no need to extend this.
%                   if mean(coverageRatios) >= maxCover
%                        validExtSubs(subItr) = 0;
%                   end
%              else
%                   validInstances = shareabilityArr <= maxShareability;
%              end
             validInstances = shareabilityArr <= maxShareability;

 %            % If there are full instances that do not extension (cover
 %            % enough of RF), we delete them.
 %            if mean(coverageRatios) < minRFCoverage
 %                validSubs(subItr) = 0;
 %            end
             
             %% If the coverage is too small, we don't consider this sub.
             if nnz(validInstances) == 0
                 validSubs(subItr) = 0;
             end
             
              %% Update instances to keep only valid ones (which cover most of RF).
              sub.instanceCenterIdx = sub.instanceCenterIdx(validInstances,:);
              sub.instanceChildren = sub.instanceChildren(validInstances,:);
              sub.instanceSigns = sub.instanceSigns(validInstances,:);
        end
        
        % We compress the object graph using the children, and the
        % edges they are involved. 
        subs(subItr) = sub;
        
        %% Assign the score of the sub, as well as its normalized mdl score if applicable.
        subs(subItr).mdlScore = double(subScore);
    end
end