%> Name: getPartStats
%>
%> Description: Given a set of parts, this function calculates certain
%> properties of the parts in the list. Fscore, precision, recall, assigned
%> class, shareability, coverage properties are calculated.
%> 
%> @param bestSubs Initial set of substructures.
%> @param categoryArrIdx Class ids of each image.
%> @param imageIdx The image indices of immediate children.
%> @param remainingChildren Remaining immediate children.
%> @param allLeafNodes Leaf nodes represented by each immediate child.
%> @param level1Coords Downsampled coords for layer 1.
%>
%> @retval partStas A struct that includes statistics for every part.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 22.09.2016
function [ partStats, allFeatures, assignedClassArr ] = calculatePartStats(bestSubs, categoryArrIdx, imageIdx, allLeafNodes, numberOfRemainingLeafNodes)
     %% Pre-processing.
    % Initialize data structures.
    categoryArrIdx = categoryArrIdx';
    numberOfBestSubs = numel(bestSubs);
    numberOfImages = numel(categoryArrIdx);
    numberOfCategories = numel(unique(categoryArrIdx));
    
    % Put the instances ([labelId, imageId] triple) into an array.
    allInstances = cell(numel(bestSubs),1);
    centerIdxArr = {bestSubs.instanceCenterIdx};
    parfor bestSubItr = 1:numel(bestSubs)
         relevantCenters = centerIdxArr{bestSubItr};
         allInstances{bestSubItr} = cat(2, imageIdx(relevantCenters),...
           bestSubItr * ones(numel(relevantCenters), 1));
    end
    allInstances = cat(1, allInstances{:});
    allInstances = sortrows(allInstances);
    
    % Create the feature vectors to be classified.
    [imageIds, orderIdx, ~] = unique(allInstances(:,1), 'R2012a');
    allFeatures = sparse(numel(imageIds), numberOfBestSubs) > 0;
    orderIdx = cat(1, orderIdx, (size(allInstances,1)+1));
    
    % Learn train/validation features, and assign labels.
    for imageItr = 1:numel(imageIds)
        instanceImageIdx = allInstances(orderIdx(imageItr):(orderIdx(imageItr+1)-1), 2);
        allFeatures(imageItr, fastsortedunique(instanceImageIdx)) = 1;
    end
    newFeatures = sparse(numberOfImages, numberOfBestSubs) > 0;
    newFeatures(imageIds, :) = allFeatures;
    allFeatures = newFeatures;
    clear newFeatures;
    
    % Assign train/test features and normalize them (converting them to
    % unit vectors).
    allLabels = categoryArrIdx;
    
    %% Calculate discriminative features.
    fscoreArr = zeros(numberOfBestSubs,1, 'single');
    precisionArr = zeros(numberOfBestSubs,1, 'single');
    recallArr = zeros(numberOfBestSubs,1, 'single');
    assignedClassArr = zeros(numberOfBestSubs,1);
    parfor bestSubItr = 1:numberOfBestSubs
         features = full(allFeatures(:,bestSubItr));
         assignedLabels = allLabels(features);
         assignedClass = mode(double(assignedLabels));
         assignedClassArr(bestSubItr) = assignedClass;
         precision = nnz(assignedLabels == assignedClass) / numel(assignedLabels);
         recall = nnz(assignedLabels == assignedClass) / nnz(allLabels == assignedClass);
         fscore = 2 * precision * recall / (precision+recall);
         fscoreArr(bestSubItr) = single(fscore);
         precisionArr(bestSubItr) = single(precision);
         recallArr(bestSubItr) = single(recall);
    end
    
    % Add them to the return structure.
    partStats.fscoreArr = fscoreArr;
    partStats.precisionArr = precisionArr;
    partStats.recallArr = recallArr;
    
    %% Find out how unique each part is in each image.
    uniquenessArr = ones(numberOfBestSubs,1, 'single');
    orderedInstances = sortrows(allInstances(:,[2,1])); 
    [~, IA, ~] = unique(orderedInstances(:,1), 'legacy');
    IA = [0; IA];
    parfor bestSubItr = 1:numberOfBestSubs
         instanceImageIds = allInstances((IA(bestSubItr)+1):IA(bestSubItr+1),1);
         uniquenessArr(bestSubItr) = numel(fastsortedunique(sort(instanceImageIds))) / numel(instanceImageIds);
    end
    partStats.uniquenessArr = uniquenessArr;
    
    %% Find out shareability.
    uniqueInstances = unique(allInstances, 'rows');
    % Convert to partId, imageId order.
    uniqueInstances = uniqueInstances(:,[2,1]);
    uniqueInstances = sortrows(uniqueInstances);
    [~, IA, ~] = unique(uniqueInstances(:,1), 'legacy');
    IA = [0; IA];
    shareabilityArr = zeros(numberOfBestSubs,1, 'single');
    sampleShareabilityArr = zeros(numberOfBestSubs,1, 'single');
    parfor bestSubItr = 1:numberOfBestSubs
         relevantInstances = uniqueInstances((IA(bestSubItr)+1):IA(bestSubItr+1),:);
         shareabilityArr(bestSubItr) = single(numel(fastsortedunique(sort(categoryArrIdx(relevantInstances(:,2))))) / numberOfCategories);
         sampleShareabilityArr(bestSubItr) = single(numel(relevantInstances(:,1)) / numberOfImages);
    end
    partStats.shareabilityArr = shareabilityArr;
    partStats.sampleShareabilityArr = sampleShareabilityArr;
    
    %% Calculate reconstructive features.
   % Allocate space for log likelihood results.
   subCoveredNodes = cell(numberOfBestSubs,1);
   
   % Go over all possible part-subpart pairs, and calculate probabilities.
   parfor subItr = 1:numberOfBestSubs
       instanceChildren = bestSubs(subItr).instanceChildren;
       instanceChildren = instanceChildren(:);
       
        % Save child probabilities.
        allNodes = fastsortedunique(sort(instanceChildren));
        coveredNodes = fastsortedunique(sort(cat(2, allLeafNodes{allNodes}))); %#ok<PFBNS>
        subCoveredNodes{subItr} = coveredNodes;
   end
   partStats.coverageArr = single(cellfun(@(x) numel(x)/numberOfRemainingLeafNodes, subCoveredNodes));
end

