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
    categoryArrIdx, allFeatures, assignedClassArr, fscoreArr)

    if numel(bestSubs) < numberOfFinalSubs
       numberOfFinalSubs = numel(bestSubs); 
    end
    
    % Initialize data structures.
    categoryArrIdx = categoryArrIdx';
    allLabels = categoryArrIdx;
    
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
%     opt = statset('display','iter');
%     cvp = cvpartition(allLabels','k',5);
%     [fs, history] = sequentialfs(@classf,allFeatures,allLabels','cv',cvp,'options',opt);
    
    discriminativeSubs = mrmr_miq_d(allFeatures, allLabels, numberOfFinalSubs);
    discriminativeSubs = sort(discriminativeSubs)';
end

function err = classf(xtrain,ytrain,xtest,ytest)
         yfit = classify(xtest,xtrain,ytrain,'diagQuadratic');
         err = sum(~strcmp(ytest,yfit));
end