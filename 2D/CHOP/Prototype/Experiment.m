function [] = Experiment(datasetName, fileType)
    generateAutoFilters(datasetName, fileType);
    runVocabularyLearning(datasetName, fileType, '');
    MarkCategoryLabels(datasetName);
    runTestInference(datasetName, fileType);
    EvaluateCategorization(datasetName);
end