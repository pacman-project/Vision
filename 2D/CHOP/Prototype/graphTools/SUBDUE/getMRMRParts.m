%> Name: getMRMRParts
%>
%> Description: Given a set of parts in bestSubs, this function greedily
%> selects a set of parts that minimize the likelihood of the data. The data
%> is grouped into overlapping receptive fields, and the reduction in the
%> cost is associated with increasing likelihood of the underlying data. Two
%> factors are contributing towards the data likelihood description, namely
%> node label and position prediction. 
%> 
%> @param bestSubs Initial set of substructures.
%> @param realNodeLabels The real labels of underlying data. 
%> @param realEdgeLabels The real labels of the edges that encode spatial
%> distributions in the bottom level. 
%> @param allEdges All edges encoded in the first level, with each cell
%> corresponding to a separate node's edges. 
%> @param allEdgeProbs Probabilities associated with edges.
%> @param numberOfFinalSubs Selection will stop if the number of selected
%> subs exceeds numberOfFinalSubs. 
%> @param stoppingCoverage The minimum coverage that is required to stop
%> selection.
%> @param uniqueChildren The ids of the nodes(data) to be covered.
%>
%> @retval discriminativeSubs Ids of final subs.
%> @retval overallCoverage The coverage on the data.
%> @retval dataLikelihood Data likelihood given the selected parts.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 08.10.2015

function [discriminativeSubs, fscoreSubs] = getMRMRParts(bestSubs, numberOfFinalSubs, ...
    categoryArrIdx, imageIdx)

    if numel(bestSubs) < numberOfFinalSubs
       numberOfFinalSubs = numel(bestSubs); 
    end
    
    % Initialize data structures.
    categoryArrIdx = categoryArrIdx';
    numberOfBestSubs = numel(bestSubs);
    
    % Put the instances ([labelId, imageId] triple) into an array.
    allInstances = cell(numel(bestSubs),1);
    for bestSubItr = 1:numel(bestSubs)
       allInstances{bestSubItr} = [imageIdx(bestSubs(bestSubItr).instanceCenterIdx),...
           repmat(bestSubItr, numel(bestSubs(bestSubItr).instanceCenterIdx), 1)];
    end
    allInstances = cat(1, allInstances{:});
    allInstances = sortrows(allInstances);
    
    % Create the feature vectors to be classified.
    allFeatures = zeros(max(imageIdx), numberOfBestSubs);
    imageIds = unique(imageIdx);
    
    % Learn train/validation features, and assign labels.
    for imageItr = imageIds'
        instanceImageIdx = allInstances(allInstances(:,1) == imageItr, 2);
        featureArr = hist(instanceImageIdx, 1:numberOfBestSubs);
        allFeatures(imageItr, :) = featureArr;
    end 
    
    % Assign train/test features and normalize them (converting them to
    % unit vectors).
    validRows = sum(allFeatures,2) ~= 0;
    allFeatures(validRows,:) = normr(allFeatures(validRows,:));
    allFeatures = uint8(round(allFeatures * 255));
    allLabels = categoryArrIdx;
    
    %% Calculate f-scores of the parts so we can verify if we've selected good parts.
    fscoreArr = zeros(numel(bestSubs),1);
    assignedClassArr = zeros(numel(bestSubs),1);
    for bestSubItr = 1:numel(bestSubs)
         features = allFeatures(:,bestSubItr) > 0;
         assignedLabels = allLabels(features);
         assignedClass = mode(assignedLabels);
         assignedClassArr(bestSubItr) = assignedClass;
         precision = nnz(assignedLabels == assignedClass) / numel(assignedLabels);
         recall = nnz(assignedLabels == assignedClass) / nnz(allLabels == assignedClass);
         fscore = 2 * precision * recall / (precision+recall);
         fscoreArr(bestSubItr) = fscore;
    end
    
    %% Select parts based on f-scores.
    [~, fscoreArrSortIdx] = sort(fscoreArr, 'descend');
    sortedAssignedClassArr = assignedClassArr(fscoreArrSortIdx);
    categories = unique(allLabels);
    fscoreSubs = [];
    for categoryItr = categories
       % Get top N.
       idx = find(sortedAssignedClassArr == categoryItr, ceil(numberOfFinalSubs/numel(categories)), 'first');
       fscoreSubs = [fscoreSubs; fscoreArrSortIdx(idx)]; %#ok<AGROW>
    end
    fscoreSubs = sort(fscoreSubs);
    
    %% Here, we apply MR-MR based feature selection.
    discriminativeSubs = mrmr_miq_d(allFeatures, allLabels, numberOfFinalSubs);
    discriminativeSubs = sort(discriminativeSubs)';
    
%     % Get train and validation accuracy for final evaluation.
%     % Finally, we have to select a subset of these features since we would
%     % like to maximize the performance with the least number of possible
%     % subs.
%     maxAcc = 0;
%     maxPrec = 0;
%     maxAccItr = numel(discriminativeSubs);
%     if valItr ~= -1
%         accArr = zeros(numel(20:20:numel(discriminativeSubs)), 1);
%         for subItr = 20:20:numel(discriminativeSubs)
%             [tempAccuracy, tempPrecision] = calculateCategorizationAccuracy(bestSubs(discriminativeSubs(1:subItr)), ...
%                 categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision, 1, true);
%             accArr(round(subItr/20)) = tempAccuracy;
%             if tempAccuracy > maxAcc
%                 maxAccItr = subItr;
%                 maxPrec = tempPrecision;
%                 maxAcc = tempAccuracy;
%             end
%         end
%     end
%    
%    figure, plot(20:20:numel(discriminativeSubs), accArr), hold on;
%    axis([20 numel(discriminativeSubs) 0 1]);
%    hold off;
%     discriminativeSubs = discriminativeSubs(1:maxAccItr);
%     [trueAccuracy, truePrecision] = calculateCategorizationAccuracy(bestSubs(discriminativeSubs(1:maxAccItr)), ...
%        categoryArrIdx, imageIdx, validationIdx, valItr, midThr, singlePrecision, 1, true);
%     display(['[SUBDUE] [Val ' num2str(valItr) '] Val precision after discriminative part selection : %' ...
%         num2str(100 * truePrecision) ...
%         ', having ' num2str(numel(discriminativeSubs)) ' subs.']);  
%     
%     % Return correct indices now.
%     discriminativeSubs = validSubIdx(discriminativeSubs);
end