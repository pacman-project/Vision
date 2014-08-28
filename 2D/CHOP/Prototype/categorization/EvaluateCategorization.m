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
    options = SetParameters(datasetName, 'train');
    outputFolder = options.testInferenceFolder;
    fileNames = fuf([outputFolder '/*.mat'], 1, 'detail');
    gtArr = NaN(size(fileNames,1), 1);
    detectionArr = NaN(size(fileNames,1), 1);
    estimatedCategoryLabel = NaN;
    w = warning('off', 'all');
    for fileItr = 1:numel(fileNames)
        
       load(fileNames{fileItr}, 'categoryLabel', 'estimatedCategoryLabel');
       % If this file has not been processed yet, move on.
       if isnan(estimatedCategoryLabel)
           continue;
       else
           1;
       end
       gtArr(fileItr) = categoryLabel;
       detectionArr(fileItr) = estimatedCategoryLabel;
       estimatedCategoryLabel = NaN;
    end
    warning(w);
    
    % Find number of processed samples.
    processedSamples = ~isnan(detectionArr);
    numberOfProcessedSamples = numel(find(processedSamples));
    detectionArr = detectionArr(processedSamples);
    gtArr = gtArr(processedSamples);
    
    % Estimate performance, and find confusion matrix.
    accuracy = numel(find(detectionArr == gtArr)) / numberOfProcessedSamples;
    confMat = confusionmat(gtArr, detectionArr);
    save([options.outputFolder '/perf.mat'], 'accuracy', 'confMat');
    display(['Accuracy: ' num2str(accuracy)]);
    display('Confusion matrix:');
    display(mat2str(confMat));
    confMatStr = mat2str(confMat);
    fid = fopen([options.outputFolder '/confMat.txt'], 'w');
    fprintf(fid, '%s', confMatStr);
    fclose(fid);
end

