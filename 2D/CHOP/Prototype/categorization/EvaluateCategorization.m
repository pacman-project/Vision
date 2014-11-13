%> Name: EvaluateCategorization
%>
%> Description: Evaluates the categorization performance of the dataset.
%>
%> @param datasetName The name of the dataset.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 04.07.2014
function [ ] = EvaluateCategorization( datasetName, perfType, minLevels, maxLevels )
    % Read vocabulary.
    options = SetParameters(datasetName, 'train');
    load([options.outputFolder '/vb.mat']);
    outputFolder = options.testInferenceFolder;
    vocabulary = vocabulary;
    if strcmp(perfType, 'test') 
        fileNames = fuf([outputFolder '/*.mat'], 1, 'detail');
        gtArr = NaN(size(fileNames,1), 1);
        detectionArr = NaN(size(fileNames,1), 1);
        w = warning('off', 'all');
        for fileItr = 1:numel(fileNames)
           estimatedCategoryLabel = NaN;
           load(fileNames{fileItr});
           % If this file has not been processed yet, move on.
           if isnan(estimatedCategoryLabel)
               estimatedCategoryLabel = getCategoryLabel(vocabulary, exportArr, confidenceArr, minLevels, maxLevels);
           end
           gtArr(fileItr) = categoryLabel;
           detectionArr(fileItr) = estimatedCategoryLabel;
        end
        warning(w);
    else
        confidenceArr = [];
        load([options.outputFolder '/export.mat']);
        numberOfImages = max(exportArr(:,5));
        gtArr = categoryArrIdx;
        detectionArr = NaN(numberOfImages,1);
        for imageItr = 1:numberOfImages
            exportArrImg = exportArr(exportArr(:,5) == imageItr,:);
            detectionArr(imageItr) = getCategoryLabel(vocabulary, exportArrImg, confidenceArr, minLevels, maxLevels);
        end
    end
    
    if ~isempty(strfind(datasetName, 'MNIST'))
        categoryNames = cellfun(@(x) str2double(x), categoryNames);
        customOrder = zeros(size(categoryNames,1),1);
        customOrder(categoryNames+1) = 1:size(categoryNames,1);
    else
        customOrder = [];
    end
    
    % Find number of processed samples.
    processedSamples = ~isnan(detectionArr);
    numberOfProcessedSamples = numel(find(processedSamples));
    detectionArr = detectionArr(processedSamples);
    gtArr = gtArr(processedSamples);
    
    % Estimate performance, and find confusion matrix.
    accuracy = numel(find(detectionArr == gtArr)) / numberOfProcessedSamples;
    confMat = confusionmat(gtArr, detectionArr);
    if ~isempty(customOrder)
        confMat= confMat(customOrder, customOrder);
    end
    
    save([options.outputFolder '/perf.mat'], 'accuracy', 'confMat');
    display(['Accuracy: ' num2str(accuracy)]);
    display('Confusion matrix:');
    display(mat2str(confMat));
    confMatStr = mat2str(confMat);
    fid = fopen([options.outputFolder '/confMat.txt'], 'w');
    fprintf(fid, '%s', confMatStr);
    fclose(fid);
end

