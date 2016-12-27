%> Name: MarkCategoryLabels
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
function [ ] = MarkCategoryLabels( datasetName )
    %#ok<*NODEF>
    options = SetParameters(datasetName, 'train');
    % Read the vocabulary and the exported realizations. 
    load([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary', 'categoryNames');
    load([pwd '/output/' datasetName '/export.mat'], 'categoryArrIdx');
    load([pwd '/output/' datasetName '/export.mat'], 'exportArr');
%    load([pwd '/output/' datasetName '/export.mat'], 'categoryArrIdx', 'trainingFileNames');
%    numberOfImages = numel(categoryArrIdx);
%     
%     % Obtain categoryArrIdx and exportArr from the data.
%      allExportArr = cell(numberOfImages,1);
%     for imageItr = 1:numberOfImages
%          [~, fileName, ~] = fileparts(trainingFileNames{imageItr});
%          load([options.testInferenceFolder '/' categoryNames{categoryArrIdx(imageItr)} '_' fileName '_test.mat'], 'exportArr');
%          exportArr(:,5) = imageItr;
%          allExportArr{imageItr} = exportArr;
%     end
%     exportArr = cat(1, allExportArr{:});
    
    % We go through each layer of the vocabulary, and update categoryArr of
    % every node in the vocabulary with the probabilities that this node
    % belongs to each category. The probabilities of each category for a
    % single node sum up to 1. categoryArr of a node is 1xN array, where N
    % is the number of categories.
    for levelItr = 1:numel(vocabulary)
         
        vocabLevel = vocabulary{levelItr};
        
        if isempty(vocabLevel)
             break;
        end
        
        % Put exported realizations in different cells in a cell array 
        % (for fast parallel processing).
        realizations = exportArr(exportArr(:,4) == levelItr, :); 
        numberOfParts = numel(vocabLevel);
        numberOfCategories = numel(categoryNames);
        realArr = cell(numberOfParts,1);
        accHist = ones(1, numberOfCategories);
        
        for partItr = 1:numberOfParts
            realArr{partItr} = categoryArrIdx(unique(realizations(realizations(:,1) == partItr, 5)));
            accHist = accHist + hist(realArr{partItr}, 1:numberOfCategories);
        end
        normConstants = 1./accHist;
        normConstants = normConstants / sum(normConstants);
        
%         % Find normalization constants for each category.
%         realCategoryArr = categoryArrIdx(realizations(:, 5));
%         normConstants = hist(realCategoryArr, 1:numberOfCategories);
%         normConstants = normConstants/sum(normConstants);
%         normConstants = 1./normConstants;
        
        % Process each node in the vocabulary.
        probArr = cell(numberOfParts,1);
        for partItr = 1:numel(vocabLevel)
            if ~isempty(realArr{partItr})
                assgnArr = single((hist(realArr{partItr}, 1:numberOfCategories) / numel(realArr{partItr})));
            else
                assgnArr = zeros(1, numberOfCategories, 'single');
            end
            assgnArr = assgnArr .* normConstants;
            assgnArr(isnan(assgnArr)) = 0;
            assgnArr = assgnArr / sum(assgnArr);
            assgnArr(isnan(assgnArr)) = 0;
            probArr{partItr} = assgnArr;
        end
        
        % Save probabilities.
        [vocabLevel.categoryArr] = deal(probArr{:});
        vocabulary{levelItr} = vocabLevel; %#ok<AGROW>
    end
    
    % Save vocabulary.
    save([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary', '-append');
end

