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
    gtArr = zeros(size(fileNames,1), 1);
    detectionArr = zeros(size(fileNames,1), 1);
    for fileItr = 1:numel(fileNames)
       load(fileNames{fileItr}, 'categoryLabel', 'estimatedCategoryLabel');
       gtArr(fileItr) = categoryLabel;
       detectionArr(fileItr) = estimatedCategoryLabel;
    end
    
    % Estimate performance, and find confusion matrix.
    accuracy = numel(find(detectionArr == gtArr)) / numel(fileNames);
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

