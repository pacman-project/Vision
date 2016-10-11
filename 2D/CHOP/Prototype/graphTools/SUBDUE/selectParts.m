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
function [bestSubs, currentAccuracy] = selectParts(bestSubs, ...
    numberOfFinalSubs, categoryArrIdx, imageIdx, ...
    allSigns, allLeafNodes, level1Coords, supervisionFlag, folderName)

    % Keep a list of positive children, if required.
    posNodes = find(allSigns);

    % Warn the user. This may take a while.
    if supervisionFlag
        display('[SUBDUE] Running supervised part selection. This may take a while..');
    else
        display('[SUBDUE] Running unsupervised part selection. This may take a while..');
    end

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

    optimalCount = numberOfFinalSubs;
    currentAccuracy = -1;
    
    % Select remaining children and filter negative nodes.
    remainingChildren = fastsortedunique(sort(cat(1, allChildren{:})));
    if ~supervisionFlag
        remainingChildren = intersect(remainingChildren, posNodes);
    end
    
    % Measure number of remaining leaf nodes.
    numberOfRemainingLeafNodes = numel(fastsortedunique(sort(cat(2, allLeafNodes{remainingChildren}))));

    %% We analyse the parts from multiple viewpoints and generate a report.
    % Fscore.
    % Precision
    % Recall
    % Reconstructability
    % Coverage
    % Number of Instances
    % Number of Leaf nodes
    % Shareability
    [partStats, allFeatures, assignedClassArr] = calculatePartStats(bestSubs, categoryArrIdx, imageIdx, allLeafNodes, numberOfRemainingLeafNodes);
    
    %% Select subs.
    display('[SUBDUE] 1/2. Reconstructive part selection running...');
    [reconstructiveSubs] = getReconstructiveParts(bestSubs, optimalCount, level1Coords, remainingChildren, allLeafNodes);
    display('[SUBDUE] 2/2. Discriminative part selection running...');
    [discriminativeSubs, fscoreSubs] = getMRMRParts(bestSubs, numel(reconstructiveSubs), categoryArrIdx, allFeatures, assignedClassArr, partStats.fscoreArr);
    if supervisionFlag
        validSubs = discriminativeSubs;
    else
        validSubs = reconstructiveSubs;
    end
    
    % Visualize all part statistics.
    visualizePartStats(partStats, reconstructiveSubs, discriminativeSubs, fscoreSubs, folderName);
    
    % Update instance information.
    bestSubs = bestSubs(validSubs);
end