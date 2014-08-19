function [] = Experiment(datasetName, fileType)
    generateAutoFilters(datasetName);
    runVocabularyLearning(datasetName, fileType, '');
    MarkCategoryLabels(datasetName);
    runTestInference(datasetName, fileType);
    EvaluateCategorization(datasetName);
end