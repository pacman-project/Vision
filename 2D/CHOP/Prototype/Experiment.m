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
    diary([pwd '/logs/' datasetName '_log.txt']);
    addpath([pwd '/utilities']);
    generateAutoFilters(datasetName, fileType);
    runVocabularyLearning(datasetName, fileType, '');
    MarkCategoryLabels(datasetName);
    runTestInference(datasetName, fileType);
    EvaluateCategorization(datasetName);
    diary off;
end