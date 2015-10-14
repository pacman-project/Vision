%> Name: selectParts
%>
%> Description: Given a set of parts in bestSubs, this function 
%> calculates the optimal set of parts on the validation data, if it
%> exists. If it does not exist, training data is used. There are two 
%> evaluation mechanisms for parts: First, we evaluate them using a
%> reconstruction-based criterion. However, if supervisionFlag is true, or
%> our part selection has reduced the categorization performance on the
%> dataset, we switch to discrimination-based part selection.
%> Reconstruction-based selection: Each part is
%> evaluated based on its coverage on the data, and total matching cost of
%> the part's instances. We're trying to maximize the coverage on the
%> training set, while minimizing the matching cost of the instances. Only
%> a limited set of subs is selected. 
%> Discrimination-based selection: We categorize the images, and calculate 
%> the overall accuracy. If the overall accuracy has dropped since the last
%> level, we switch to discrimination, and set supervisionFlag as true.
%> 
%> @param bestSubs A set of substructures evaluated on the training set.
%> @param instanceChildrenDescriptors N_instance x maxNumberOfChildren
%> array that has an ordered list of children from the previous level for
%> each instance.
%> @param remainingInstanceLabels N_instance x 1 array that marks the part
%> label for each instance.
%> @param allLeafNodes N x 1 cell array that keeps what leaf nodes are covered 
%> by each part's instances.  
%> prevGraphNodeCount Number of nodes in the previous level's object
%> graphs.
%> @param stoppingCoverage Selection will stop if stoppingCoverage percent
%> of previous levels' leaf nodes are already covered.
%> @param numberOfFinalSubs Selection will stop if the number of selected
%> subs exceeds numberOfFinalSubs. 
%>
%>
%> @retval bestSubs Final set of best substructures.
%> @retval instanceChildrenDescriptors N_instance x maxNumberOfChildren
%> array that has an ordered list of children from the previous level for
%> each instance.
%> @retval remainingInstanceLabels N_instance x 1 array that marks the part
%> label for each instance.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 13.04.2015
function [bestSubs, optimalThreshold, optimalAccuracy] = selectParts(bestSubs, realNodeLabels, ...
    nodeDistanceMatrix, edgeDistanceMatrix, nodePositions, edgeCoords, ...
    singlePrecision, numberOfFinalSubs, fixedThreshold, ...
    validationFolds, validationIdx, categoryArrIdx, imageIdx, ...
    allSigns, supervisionFlag)

    % Keep a list of positive children, if required.
    posNodes = find(allSigns);

    % Warn the user. This may take a while.
    if supervisionFlag
        display('[SUBDUE] Running supervised part selection. This may take a while..');
    else
        display('[SUBDUE] Running unsupervised part selection. This may take a while..');
    end

    % Keep best subs.
    orgBestSubs = bestSubs;
    
    % Start processing with the remaining subs.
    numberOfOrgBestSubs = numel(bestSubs);
    allChildren = cell(numberOfOrgBestSubs,1);
    for bestSubItr = 1:numel(bestSubs)
        children = bestSubs(bestSubItr).instanceChildren;

        % Get unique children.
        if ~isempty(children)
              children  = children(:);
              children = fastsortedunique(sort(children));
        else 
              children = [];
        end

        % Make the sizes uniform, i.e. Nx1.
        if size(children,2) ~= 1
            children = children';
        end

        % Assign the children.
        allChildren{bestSubItr} = children;
    end

    optimalThreshold = fixedThreshold;
    optimalCount = numberOfFinalSubs;
    optimalAccuracy = -1;
    
    %% Finally, given the optimal threshold, we select the best subs based on different validation sets and aggregate them.
    aggregatedSubs = cell(validationFolds,1);
    for valItr = 1:validationFolds
        % We exclude the subs which have zero-cost matchs on this subset,
        % but not on other subsets.
        validSubIdx = ones(numel(orgBestSubs),1) > 0;

        % If there's only one fold, or we do not wish to do
        % cross-validation, we fit the training data.
        if validationFolds > 1
            for subItr = 1:numel(orgBestSubs)
                instanceMatchCosts = orgBestSubs(subItr).instanceMatchCosts;
                instanceValidationIdx = orgBestSubs(subItr).instanceValidationIdx;

                if nnz(instanceMatchCosts == 0 & instanceValidationIdx ~= valItr) == 0
                   validSubIdx(subItr) = 0; 
                end
            end
        end
        bestSubs = orgBestSubs(validSubIdx);
        numberOfBestSubs = numel(bestSubs);
        if numberOfBestSubs == 0
           continue; 
        end

        % Select remaining children and filter negative nodes.
        remainingChildren = fastsortedunique(sort(cat(1, allChildren{validSubIdx}))); %#ok<PFBNS>
        if ~supervisionFlag
            remainingChildren = intersect(remainingChildren, posNodes);
        end
        
        validSubIdx = find(validSubIdx);
        
        % Select subs.
        if supervisionFlag
            [validSubs, ~, ~] = getMRMRParts(bestSubs, optimalCount, nodeDistanceMatrix, edgeDistanceMatrix, ...
                categoryArrIdx, imageIdx, validationIdx, valItr, optimalThreshold, singlePrecision);
        else
           [validSubs, ~, ~] = getReconstructiveParts(bestSubs, realNodeLabels, ...
                nodePositions, edgeCoords, ...
                optimalCount, remainingChildren);
        end
        aggregatedSubs{valItr} = validSubIdx(validSubs);
    end
    
    % Finally, obtain a list of final subs and get their union.
    finalSubList = unique(cat(1, aggregatedSubs{:}));
    bestSubs = orgBestSubs;
    
   % Update instance information.
   bestSubs = bestSubs(finalSubList);
   % Update bestSubs instances by taking the new threshold into account.
   parfor bestSubItr = 1:numel(bestSubs)
        sub = bestSubs(bestSubItr);
        validInstances = sub.instanceMatchCosts < ((size(sub.edges,1)*2+1) * optimalThreshold + singlePrecision);
        sub.instanceCenterIdx = sub.instanceCenterIdx(validInstances,:);
        sub.instanceChildren = sub.instanceChildren(validInstances,:);
        sub.instanceMappings = sub.instanceMappings(validInstances,:);
        if ~isempty(sub.edges)
            sub.instanceEdges = sub.instanceEdges(validInstances,:);
        end
        sub.instanceSigns = sub.instanceSigns(validInstances,:);
        sub.instanceCategories = sub.instanceCategories(validInstances,:);
        sub.instanceMatchCosts = sub.instanceMatchCosts(validInstances,:);
        sub.instanceValidationIdx = sub.instanceValidationIdx(validInstances,:);
        bestSubs(bestSubItr) = sub;
   end
end