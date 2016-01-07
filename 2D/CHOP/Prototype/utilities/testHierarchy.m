%> Name: testHierarchy
%>
%> Description: This function checks the consistency of the inference
%> procedure in training and testing, and reports on the results. It's a test
%> function that should run regularly to check if the system is working
%> properly. It assumes that the training is already performed on the
%> dataset specified with 'datasetName'.
%>
%> @param datasetName Name of the dataset to check.
%>
%> @retval flag 0 if all ok, -1 if detections are wrong, -2 if precise
%> positions are wrong, -3 if activations are wrong.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.11.2015
function [ flag ] = testHierarchy( datasetName, fileExt )
     % First, we perform inference for the images in training set, if necessary. 
     inputFolder = [pwd '/input/' datasetName];
     vocabFolder  = [inputFolder '/vocab'];
     testFolder = [inputFolder '/test'];
     testTempFolder = [inputFolder '/testTemp'];
     inferenceFolder = [pwd '/output/' datasetName '/test/inference']; 
     tempInferenceFolder = [pwd '/output/' datasetName '/test/inferenceTemp'];
     flag = 0;
     
     % Replace test folder with the vocabfolder. 
      movefile(testFolder, testTempFolder);
      movefile(vocabFolder, testFolder);
      movefile(inferenceFolder, tempInferenceFolder);
     
     % Run inference.
      runTestInference(datasetName, fileExt);
     
     % One by one, we check if each image's parse tree matches in training
     % inference and test inference.
     load([pwd '/output/' datasetName '/export.mat']);
     allExportArr = exportArr;
     allPrecisePositions = precisePositions;
     allActivationArr = activationArr;
     for imgItr = 1:numel(trainingFileNames)
          [~, fileName, ~] = fileparts(trainingFileNames{imgItr});
          categoryStr = categoryNames{categoryArrIdx(imgItr)};
          fileStr = [categoryStr '_' fileName];
          load([inferenceFolder '/' fileStr '_test.mat']);
          
          % Compare exported activations, their precise positions, and
          % activation values.
          idx = allExportArr(:,5) == imgItr;
          exportArrTrain = allExportArr(idx,:);
          exportArrTrain(:,5) = 1;
          precisePositionsTrain = allPrecisePositions(idx,:);
          activationArrTrain = allActivationArr(idx,:);
          
          if ~isequal(exportArrTrain, exportArr)
               flag = -1;
               break;
          end
          if ~isequal(round(precisePositionsTrain*1000) / 1000, round(precisePositions* 1000) / 1000)
               flag = -2;
          end
          if ~isequal(activationArrTrain, activationArr)
               flag = -3;
          end
     end
     
     % Move back the folders.
      movefile(testFolder, vocabFolder);
      movefile(testTempFolder, testFolder);
      movefile(tempInferenceFolder, inferenceFolder);
end