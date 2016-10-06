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
function [ partStats ] = calculatePartStats(bestSubs, categoryArrIdx, imageIdx, allLeafNodes)
     %% Pre-processing.
    % Initialize data structures.
    categoryArrIdx = categoryArrIdx';
    numberOfBestSubs = numel(bestSubs);
    numberOfImages = numel(categoryArrIdx);
    numberOfCategories = numel(unique(categoryArrIdx));
    
    % Put the instances ([imageId, labelId] duplet) into an array.
    allInstances = cell(numberOfBestSubs,1);
    for bestSubItr = 1:numberOfBestSubs
       allInstances{bestSubItr} = [imageIdx(bestSubs(bestSubItr).instanceCenterIdx),...
           repmat(bestSubItr, numel(bestSubs(bestSubItr).instanceCenterIdx), 1)];
    end
    allInstances = cat(1, allInstances{:});
    allInstances = sortrows(allInstances);
    
    % Create the feature vectors to be classified.
    allFeatures = zeros(numberOfImages, numberOfBestSubs);
    imageIds = unique(allInstances(:,1));
    
    [~, IA, ~] = unique(allInstances(:,1), 'legacy');
    IA = [0; IA];
    % Learn train/validation features, and assign labels.
    for imageItr = 1:numel(imageIds)
        instanceImageIdx = allInstances((IA(imageItr)+1):IA(imageItr+1), 2);
        featureArr = hist(instanceImageIdx, 1:numberOfBestSubs);
        allFeatures(imageIds(imageItr), :) = featureArr;
    end 
    
    % Assign train/test features and normalize them (converting them to
    % unit vectors).
    validRows = sum(allFeatures,2) ~= 0;
    allFeatures(validRows,:) = normr(allFeatures(validRows,:));
    allFeatures = uint8(round(allFeatures * 255));
    allLabels = categoryArrIdx;
    
    %% Calculate discriminative features.
    fscoreArr = zeros(numberOfBestSubs,1, 'single');
    precisionArr = zeros(numberOfBestSubs,1, 'single');
    recallArr = zeros(numberOfBestSubs,1, 'single');
    assignedClassArr = zeros(numberOfBestSubs,1);
    for bestSubItr = 1:numberOfBestSubs
         features = allFeatures(:,bestSubItr) > 0;
         assignedLabels = allLabels(features);
         assignedClass = mode(assignedLabels);
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
    for bestSubItr = 1:numberOfBestSubs
         instanceImageIds = allInstances((IA(bestSubItr)+1):IA(bestSubItr+1),1);
         uniquenessArr(bestSubItr) = numel(unique(instanceImageIds)) / numel(instanceImageIds);
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
    for bestSubItr = 1:numberOfBestSubs
         relevantInstances = uniqueInstances((IA(bestSubItr)+1):IA(bestSubItr+1),:);
         shareabilityArr(bestSubItr) = single(numel(unique(categoryArrIdx(relevantInstances(:,2)))) / numberOfCategories);
         sampleShareabilityArr(bestSubItr) = single(numel(relevantInstances(:,1)) / numberOfImages);
    end
    partStats.shareabilityArr = shareabilityArr;
    partStats.sampleShareabilityArr = sampleShareabilityArr;
    
    %% Calculate reconstructive features.
   % Allocate space for log likelihood results.
   subCoveredNodes = cell(numberOfBestSubs,1);
   
   % Go over all possible part-subpart pairs, and calculate probabilities.
   for subItr = 1:numberOfBestSubs
       instanceChildren = bestSubs(subItr).instanceChildren;
       
        % Save child probabilities.
        allNodes = fastsortedunique(sort(instanceChildren));
        allNodes = (allNodes(:))';
        coveredNodes = fastsortedunique(sort(cat(2, allLeafNodes{allNodes})));
        subCoveredNodes{subItr} = coveredNodes;
   end
   partStats.coverageArr = single(cellfun(@(x) numel(x), subCoveredNodes));
end

