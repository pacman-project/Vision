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
    load([pwd '/output/' datasetName '/export.mat'], 'categoryArrIdx');
    load([pwd '/output/' datasetName '/trainingInference.mat'], 'exportArr', 'activationArr');

%     load([pwd '/output/' datasetName '/export.mat'], 'categoryArrIdx', 'trainingFileNames');
%     numberOfImages = numel(categoryArrIdx);
%     
%     % Obtain categoryArrIdx and exportArr from the data.
%      allExportArr = cell(numberOfImages,1);
%      allActivationArr = cell(numberOfImages,1);
%     for imageItr = 1:numberOfImages
%          [~, fileName, ~] = fileparts(trainingFileNames{imageItr});
%          load([options.testInferenceFolder '/' categoryNames{categoryArrIdx(imageItr)} '_' fileName '_test.mat'], 'exportArr', 'activationArr');
%          exportArr(:,5) = imageItr;
%          allExportArr{imageItr} = exportArr;
%          allActivationArr{imageItr} = activationArr;
%     end
%     exportArr = cat(1, allExportArr{:});
%     activationArr = cat(1, allActivationArr{:});
    numberOfImages = numel(categoryArrIdx);
    
    % Go through every image in the training set, and classify.
    confMat = zeros(numel(categoryNames));
    
    top1Mat = zeros(numberOfImages, 1) > 0;
    top3Mat = zeros(numberOfImages, 1) > 0;
    top5Mat = zeros(numberOfImages, 1) > 0;
    for imageItr = 1:numberOfImages
         % Obtain top-level realizations for this image.
         idx = exportArr(:,5) == imageItr;
         realizations = exportArr(idx,:);
         activations = activationArr(idx);
         maxLevel = max(realizations(:,4));
         idx2=realizations(:,4) == maxLevel;
         relevantRealizations = realizations(idx2,:);
         relevantActivations = activations(idx2);
         
         % Order realizations wrt activations.
         [relevantActivations, sortOrder] = sort(relevantActivations, 'descend');
         relevantRealizations = relevantRealizations(sortOrder, :);
         
         % Obtain contributing subs.
         [contributingSubs, IC, ~] = unique(relevantRealizations(:,1), 'stable');
%          if numel(contributingSubs) > 5
%               contributingSubs = contributingSubs(1:5);
%          end
         contributingSubActivations = relevantActivations(IC);
         combinedArr = cat(1, vocabulary{maxLevel}(contributingSubs).categoryArr);
         if size(combinedArr,1) > 1
              for combItr = 1:size(combinedArr,1)
                   combinedArr(combItr, :) = combinedArr(combItr, :) * (1/-contributingSubActivations(combItr));
              end
              combinedArr = mean(combinedArr, 1);
         end
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
    top1Accuracy = nnz(top1Mat) / numel(top1Mat);
    top3Accuracy = nnz(top3Mat) / numel(top3Mat);
    top5Accuracy = nnz(top5Mat) / numel(top5Mat);
    save([pwd '/output/' datasetName '/trainingClassification.mat'], 'confMat', 'top1Accuracy', 'top3Accuracy', 'top5Accuracy');
    display(['Top 1 classification performance: %' num2str(100 * top1Accuracy) '.']);
    display(['Top 3 classification performance: %' num2str(100 * top3Accuracy) '.']);
    display(['Top 5 classification performance: %' num2str(100 * top5Accuracy) '.']);
end

