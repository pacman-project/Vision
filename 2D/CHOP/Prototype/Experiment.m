%> Name: Experiment
%>
%> Description:  Given the name and file type (e.g. '.png') of the data,
%this function runs vocabulary learning and test inference codes.
%>
%> @param datasetName Name of the dataset to work on. 
%> @param datasetName File type of the dataset to work on (e.g. '.png'). 
%>
%> @retval none
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 08.09.2014
function [] = Experiment(datasetName, fileType)
    if exist([pwd '/logs/' datasetName '_log.txt'], 'file')
        delete([pwd '/logs/' datasetName '_log.txt']);
    end
    diary([pwd '/logs/' datasetName '_log.txt']);
    addpath([pwd '/utilities']);
    runVocabularyLearning(datasetName, fileType, '.png');
    MarkCategoryLabels(datasetName);
    runTrainingInference( datasetName, fileType );
    MarkCategoryLabels(datasetName);
    ClassifyTrainingImages(datasetName);
    runTestInference(datasetName, fileType);
    baselineClassification(datasetName);
    
    % Run node analysis.
    try
         runDiscriminativeAnalysis(datasetName);
         plotDiscriminativeAnalysis(datasetName);
    catch
         display('Discriminative analysis failed');
    end
    diary off;
end