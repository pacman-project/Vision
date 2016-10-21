%> Name: ClassifyTrainingImages
%>
%> Description: This function classifies training images based on CHOP's
%> output. It's based on counting by each member's contribution.
%>
%> @param datasetName The name of the dataset.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 06.10.2016
function [confMat] = ClassifyTrainingImages( datasetName )
    %#ok<*NODEF>
%    options = SetParameters(datasetName, 'train');
    % Read the vocabulary and the exported realizations. 
    load([pwd '/output/' datasetName '/vb.mat'], 'vocabulary', 'categoryNames', 'options');
    load([pwd '/output/' datasetName '/export.mat'], 'categoryArrIdx', 'exportArr');
    numberOfImages = numel(categoryArrIdx);
%     load([pwd '/output/' datasetName '/export.mat'], 'categoryArrIdx', 'trainingFileNames');
%     numberOfImages = numel(categoryArrIdx);
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
    
    % Go through every image in the training set, and classify.
    confMat = zeros(numel(categoryNames));
    
    top1Mat = zeros(numberOfImages, 1) > 0;
    top3Mat = zeros(numberOfImages, 1) > 0;
    top5Mat = zeros(numberOfImages, 1) > 0;
    for imageItr = 1:numberOfImages
         % Obtain top-level realizations for this image.
         realizations = exportArr(exportArr(:,5) == imageItr,:);
         maxLevel = max(realizations(:,4));
         relevantRealizations = realizations(realizations(:,4) == maxLevel,:);
         
         % Obtain contributing subs.
         contributingSubs = unique(relevantRealizations(:,1));
         combinedArr = cat(1, vocabulary{maxLevel}(contributingSubs).categoryArr);
         combinedArr = mean(combinedArr, 1);
         [maxVal, predictedCategory] = max(combinedArr);
         
         % In case of more than 1 values, move on.
         if nnz(combinedArr == maxVal) > 1
              predictedCategory = datasample(find(combinedArr == maxVal), 1);
         end
         realCategory = categoryArrIdx(imageItr);
         confMat(realCategory, predictedCategory) = confMat(realCategory, predictedCategory) + 1;
         
         [~, categoryOrder] = sort(combinedArr, 'descend');
         %% Measure top 1, 3 and 5 classification rates.
         top1Mat(imageItr) =ismember(realCategory, categoryOrder(1));
         top3Mat(imageItr) =ismember(realCategory, categoryOrder(1:3));
         top5Mat(imageItr) =ismember(realCategory, categoryOrder(1:5));
    end
    display(['Top 1 classification performance: %' num2str(100 * nnz(top1Mat) / numel(top1Mat)) '.']);
    display(['Top 3 classification performance: %' num2str(100 * nnz(top3Mat) / numel(top3Mat)) '.']);
    display(['Top 5 classification performance: %' num2str(100 * nnz(top5Mat) / numel(top5Mat)) '.']);
end

