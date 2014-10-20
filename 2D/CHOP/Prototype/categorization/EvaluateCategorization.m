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
function [ ] = EvaluateCategorization( datasetName )
    % Read vocabulary.
    options = SetParameters(datasetName, 'train');
    load([options.outputFolder '/vb.mat']);
    outputFolder = options.testInferenceFolder;
    fileNames = fuf([outputFolder '/*.mat'], 1, 'detail');
    gtArr = NaN(size(fileNames,1), 1);
    detectionArr = NaN(size(fileNames,1), 1);
    w = warning('off', 'all');
    for fileItr = 1:numel(fileNames)
       estimatedCategoryLabel = NaN;
       load(fileNames{fileItr});
       % If this file has not been processed yet, move on.
       if isnan(estimatedCategoryLabel)
           estimatedCategoryLabel = getCategoryLabel(vocabulary, exportArr);
       end
       if estimatedCategoryLabel ~= categoryLabel
           1;
       end
       
       gtArr(fileItr) = categoryLabel;
       detectionArr(fileItr) = estimatedCategoryLabel;
    end
    warning(w);
    
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

