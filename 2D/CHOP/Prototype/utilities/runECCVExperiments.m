%% This function runs ECCV 2014 Experiments on a number of datasets. 
% Over each set, a vocabulary is learnt, and its vital output data
% structures are saved. 
function [ ] = runECCVExperiments( expType )
    ext = '.png';
    %% Get available sets and run CHOP on each.
    setArr = [1 2 3 4 5 10 15 20 30];
    for setItr = 1:numel(setArr)
        datasetName = ['ECCV2014/' expType '/' num2str(setArr(setItr))];
        runVocabularyLearning(datasetName, ext, '');
        testInferenceTime = runTestInference(datasetName, ext);
        save([pwd '/output/ECCV2014/' expType '/' num2str(setArr(setItr)) '/tetime.mat'], 'testInferenceTime');
    end
end

