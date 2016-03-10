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
function [subs, extendSubs, validSubs, validExtSubs] = evaluateSubs(subs, extendSubs, evalMetric, allEdges, allEdgeNodePairs, allSigns, graphSize, overlap, mdlNodeWeight, mdlEdgeWeight, isMDLExact, isSupervised, isMDLNormalized, ...
     allLeafNodes, minRFCoverage, possibleLeafNodes)

    maxLeafCounts = cellfun(@(x) numel(x), possibleLeafNodes);
    numberOfSubs = numel(subs);
    validSubs = ones(numberOfSubs,1) > 0;
    validExtSubs = ones(numberOfSubs,1) > 0;
    for subItr = 1:numberOfSubs
        % Find the weight of this node, by taking the max of the category distribution. 
        if isSupervised
            categoryArr = double(subs(subItr).instanceCategories);
            weight = nnz(categoryArr == mode(categoryArr)) / numel(categoryArr);
        else
            weight = 1;
        end
        
        % Calculate Substructure score.
        [subScore, sub, numberOfNonoverlappingInstances] = getSubScore(subs(subItr), allEdges, allEdgeNodePairs, evalMetric, ...
           allSigns, mdlNodeWeight, mdlEdgeWeight, ....
            overlap, isMDLExact);
       
       if ~isempty(extendSubs)
            extendSub = extendSubs(subItr);
       else
            extendSub = [];
       end
        
       % If this sub has multiple children, we process the instances to
       % make sure those instances which cover a large portion of their
       % receptive field survives.
        if size(subs(subItr).instanceChildren,2) > 1
             allChildren = subs(subItr).instanceChildren;
             % Here, we check if we're actually covering enough of every
             % receptive field out there.
             numberOfInstances = size(allChildren,1);
              coveredLeafNodes = cell(numberOfInstances,1);
              for childItr = 1:numberOfInstances
                   coveredLeafNodes{childItr} = int32(fastsortedunique(sort(cat(2, allLeafNodes{allChildren(childItr,:)}))));
              end

             % Find all covered leaf nodes.
%             allRFs = {allEdges(allChildren(:,1)).adjInfo};
%             allRFLeafNodes = cellfun(@(x) x(:,1:2), allRFs, 'UniformOutput', false);
%             allRFLeafNodes = cellfun(@(x) fastsortedunique(sort(cat(2, allLeafNodes{x(:)}))), allRFLeafNodes, 'UniformOutput', false)';

             % Find the intersection of two sets, to assess coverage.
             coveredLeafNodeCount = cellfun(@(x) numel(x), coveredLeafNodes);
             maxCoverLeafNodeCount = maxLeafCounts(allChildren(:,1));
             coverageRatios =  coveredLeafNodeCount ./ maxCoverLeafNodeCount;
             validInstances = coverageRatios >= minRFCoverage;
             
             % If there are full instances that do not extension (cover
             % enough of RF), we delete them.
%             validInstancesExt = coverageRatios < 1;
             if mean(coverageRatios) >= minRFCoverage
                  validExtSubs(subItr) = 0;
             end
             
             %% If the coverage is too small, we don't consider this sub.
             if nnz(validInstances) == 0
                  validSubs(subItr) = 0;
             end
             
             %% Update instances to keep only valid ones (which cover most of RF).
             sub.instanceCategories = sub.instanceCategories(validInstances,:);
             sub.instanceCenterIdx = sub.instanceCenterIdx(validInstances,:);
             sub.instanceChildren = sub.instanceChildren(validInstances,:);
             sub.instanceEdges = sub.instanceEdges(validInstances,:);
             sub.instanceExactMatchFlags = sub.instanceExactMatchFlags(validInstances,:);
             sub.instanceMappings = sub.instanceMappings(validInstances,:);
             sub.instanceMatchCosts = sub.instanceMatchCosts(validInstances,:);
             sub.instanceSigns = sub.instanceSigns(validInstances,:);
             sub.instanceValidationIdx = sub.instanceValidationIdx(validInstances,:);
%              
%              if ~isempty(extendSub) && nnz(validInstancesExt) > 0
%                   1
%              else
%                   2
%                   
%              end
        end
        
        % We compress the object graph using the children, and the
        % edges they are involved. 
        subScore = subScore * weight;
        subs(subItr) = sub;
        
        %% Assign the score of the sub, as well as its normalized mdl score if applicable.
        subs(subItr).mdlScore = double(subScore);
        if strcmp(evalMetric, 'mdl') || strcmp(evalMetric, 'likelihood')
            subs(subItr).normMdlScore = 1 - (double(subScore) / graphSize);
            
            % If the normalized MDL score is being asked for, we divide
            % subScore by the number of valid instances.
            if isMDLNormalized
                 subs(subItr).mdlScore = subs(subItr).mdlScore / numberOfNonoverlappingInstances;
            end
        end
    end
end