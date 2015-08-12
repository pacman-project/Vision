%> Name: CalculateConditionalProbabilities
%>
%> Description: this function works on the output of CHOP. It basically
%> lists the contributions of each node in the vocabulary to each category,
%> which helps with the categorization step, using the training images'
%> labels. The vocabulary is saved again to the same place after
%> processing.
%>
%> @param datasetName The name of the dataset.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 03.07.2014
function [ ] = runDiscriminativeAnalysis( datasetName )
    %#ok<*NODEF>
    options = SetParameters(datasetName, 'train');
    % Read the vocabulary and the exported realizations. 
    load([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary', 'categoryNames', 'distanceMatrices');
    load([options.currentFolder '/output/' datasetName '/export.mat'], 'exportArr', 'categoryArrIdx', 'poseArr');
    
    % We go through each layer of the vocabulary, and calculate precision,
    % recall, f-score of each part. 
    totalNumberOfParts = sum(cellfun(@(x) numel(x), vocabulary)); %#ok<USENS>
    precisionArr = nan(totalNumberOfParts, 2 + numel(categoryNames));
    recallArr = nan(totalNumberOfParts, 2 + numel(categoryNames));
    fscoreArr = nan(totalNumberOfParts, 2 + numel(categoryNames));
    fscoreInstanceArr = nan(totalNumberOfParts, 3);
    shareabilityArr = zeros(totalNumberOfParts, 3);
    mdlScoreArr = zeros(totalNumberOfParts, 4);
    partOffset = 1;
    for levelItr = 1:numel(vocabulary)
        
        vocabLevel = vocabulary{levelItr};
        vocabLevelLabels = cat(1, vocabLevel.label);
        [~, IA, ~] = unique(vocabLevelLabels, 'stable');
        
        % Put exported realizations in different cells in a cell array 
        % (for fast parallel processing).
        realizations = exportArr(exportArr(:,4) == levelItr, :); 
        
        % Based on pairwise distances, group parts that have a distance of
        % zero.
        realizations(:,1) = vocabLevelLabels(realizations(:,1));
        numberOfParts = numel(unique(vocabLevelLabels));
        
        %% Assign identifiers to all arrays.
        precisionArr(partOffset:(partOffset+numberOfParts-1),1) = levelItr;
        precisionArr(partOffset:(partOffset+numberOfParts-1),2) = 1:numberOfParts;
        recallArr(partOffset:(partOffset+numberOfParts-1),1:2) = ...
            precisionArr(partOffset:(partOffset+numberOfParts-1),1:2);
        fscoreArr(partOffset:(partOffset+numberOfParts-1),1:2) = ...
            precisionArr(partOffset:(partOffset+numberOfParts-1),1:2);
        fscoreInstanceArr(partOffset:(partOffset+numberOfParts-1),1:2) = ...
            precisionArr(partOffset:(partOffset+numberOfParts-1),1:2);
        shareabilityArr(partOffset:(partOffset+numberOfParts-1),1:2) = ...
            precisionArr(partOffset:(partOffset+numberOfParts-1),1:2);
        mdlScoreArr(partOffset:(partOffset+numberOfParts-1),1:2) = ...
            precisionArr(partOffset:(partOffset+numberOfParts-1),1:2);
        
        if levelItr>1
            mdlScores = cat(1, vocabLevel(IA).mdlScore);
            normMdlScores = cat(1, vocabLevel(IA).normMdlScore);
            mdlScoreArr(partOffset:(partOffset+numberOfParts-1),3:4) = [mdlScores, normMdlScores];
        end
        
        %% Calculate the metrics for every part and category pair. 
        for partItr = 1:numberOfParts
            sampleLabelIds = categoryArrIdx(unique(realizations(realizations(:,1) == partItr, 5)));
            for categoryItr = 1:numel(categoryNames)
                % Calculate precision, recall and F1 score.
               precision = nnz(sampleLabelIds == categoryItr)/numel(sampleLabelIds) * 100;
               precisionArr(partOffset + partItr - 1, 2 + categoryItr) = precision;
               recall = nnz(sampleLabelIds == categoryItr)/nnz(categoryArrIdx == categoryItr) * 100;
               recallArr(partOffset + partItr - 1, 2 + categoryItr) = recall;
               if precision > 0 && recall > 0
                   fscoreArr(partOffset + partItr - 1, 2 + categoryItr) = ...
                       (2 * precision * recall) / (precision + recall);
               end
            end
            
            % Finally, we calculate category shareability.
            shareabilityArr(partOffset + partItr - 1, 3) = (numel(unique(sampleLabelIds)) / numel(categoryNames)) * 100;
            fscoreInstanceArr(partOffset + partItr - 1, 3) = 2 * (1/numel(sampleLabelIds)) / (1 + 1/numel(sampleLabelIds));
        end
        
        % Save probabilities.
        partOffset = partOffset + numberOfParts;
    end
    
    shareabilityArr = shareabilityArr(1:(partOffset-1),:);
    recallArr = recallArr(1:(partOffset-1),:);
    fscoreArr = fscoreArr(1:(partOffset-1),:);
    precisionArr = precisionArr(1:(partOffset-1),:);
    fscoreInstanceArr = fscoreInstanceArr(1:(partOffset-1),:);
    mdlScoreArr = mdlScoreArr(1:(partOffset-1),:);
    
    % Save the data.
    if ~exist([options.currentFolder '/categorization/analysis/' datasetName], 'dir')
       mkdir([options.currentFolder '/categorization/analysis/' datasetName]);
    end
    save([options.currentFolder '/categorization/analysis/' datasetName '/discriminativeAnalysis.mat'], ...
        'precisionArr', 'recallArr', 'fscoreArr', 'shareabilityArr', 'mdlScoreArr', 'fscoreInstanceArr');
end

