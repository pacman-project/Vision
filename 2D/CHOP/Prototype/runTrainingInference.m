%> Name: runTestInference
%>
%> Description: Assuming the unsupervised vocabulary has been learned up to
%> a certain level, this function runs inference on test images to find
%> realizations of parts in vocabulary. This process is unsupervised, and
%> comes to an end in the last level in the vocabulary.
%>
%> @param datasetName Name of the dataset to work on. 
%> @param ext file extension.
%> 
%> @retval options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.02.2014
%> Ver 1.1 on 12.01.2014 Timing is added by Mete
function [] = runTrainingInference( datasetName, ext )
     vocabFolder = [pwd '/input/' datasetName '/vocab'];
     testFolder = [pwd '/input/' datasetName '/test'];
     
     % Swap folders.
     movefile(testFolder, [testFolder 'Temp']);
     movefile(vocabFolder, testFolder);

     runTestInference(datasetName, ext);
     testInferenceFolder = [pwd '/output/' datasetName '/test/inference'];
     
     % Create new export array out of the inferred parts.
     load([pwd '/output/' datasetName '/export.mat'], 'categoryArrIdx', 'categoryNames', 'trainingFileNames');
     numberOfImages = numel(categoryArrIdx);
    
     %  Obtain categoryArrIdx and exportArr from the data.
     allExportArr = cell(numberOfImages,1);
     allActivationArr = cell(numberOfImages,1);
    for imageItr = 1:numberOfImages
         [~, fileName, ~] = fileparts(trainingFileNames{imageItr});
         load([testInferenceFolder '/' categoryNames{categoryArrIdx(imageItr)} '_' fileName '_test.mat'], 'exportArr', 'activationArr');
         exportArr(:,5) = imageItr;
         allExportArr{imageItr} = exportArr;
         allActivationArr{imageItr} = activationArr;
    end
    exportArr = cat(1, allExportArr{:});
    activationArr = cat(1, allActivationArr{:});
    save([pwd '/output/' datasetName '/trainingInference.mat'], 'exportArr', 'activationArr');
    
     % Move output to training folders.
     trainingInferenceFolder = [pwd '/output/' datasetName '/training'];
     if exist(trainingInferenceFolder, 'dir')
          rmdir(trainingInferenceFolder, 's');
     end
     mkdir(trainingInferenceFolder);
     movefile(testInferenceFolder, trainingInferenceFolder);
     
     % Swap folders back.
     movefile(testFolder, vocabFolder);
     movefile([testFolder 'Temp'], testFolder);
end