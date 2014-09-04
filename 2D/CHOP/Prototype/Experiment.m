function [] = Experiment(datasetName, fileType)
    diary([pwd '/logs/' datasetName '_log.txt']);
    addpath([pwd '/utilities']);
%    generateAutoFilters(datasetName, fileType);
    runVocabularyLearning(datasetName, fileType, '');
%    MarkCategoryLabels(datasetName);
%    runTestInference(datasetName, fileType);
%    EvaluateCategorization(datasetName);
    diary off;
end