function [] = Experiment(datasetName, fileType)
    addpath([pwd '/utilities']);
    generateAutoFilters(datasetName, fileType);
    runVocabularyLearning(datasetName, fileType, '');
    MarkCategoryLabels(datasetName);
    runTestInference(datasetName, fileType);
    EvaluateCategorization(datasetName);
end