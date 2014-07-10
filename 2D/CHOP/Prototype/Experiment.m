function [] = Experiment(datasetName, fileType)
    runVocabularyLearning(datasetName, fileType, '');
    MarkCategoryLabels(datasetName);
    runTestInference(datasetName, fileType);
    EvaluateCategorization(datasetName);
end