%% This function collects statistics generated by CHOP on the set of experiments.
function [ results ] = collectECCVStatistics( expType )
    %% Get available sets and run CHOP on each.
    setArr = [1 2 3 4 5 10 15 20 30];
    topXMdl = 10;
    numberOfSets = numel(setArr);
    
    results.vocabSizeArr = zeros(1, numberOfSets);
    results.avgMdlScoreArr = zeros(1, numberOfSets);
    results.mdlScoreArr = cell(1, numberOfSets);
    results.trainTimeArr = zeros(1, numberOfSets);
    results.testTimeArr = zeros(1, numberOfSets);
    results.avgShareabilityArr = zeros(1, numberOfSets);
    results.shareabilityArr = cell(1, numberOfSets);
    
    for setItr = 1:numberOfSets
        outputFolder = [pwd '/output/ECCV2014/' expType '/' num2str(setArr(setItr))];
        
        %% Get vocabulary related features.
        load([outputFolder '/vb.mat'], 'vocabulary');
        
        % size
        vocabSizes = cellfun(@(x) numel(x), vocabulary);
        results.vocabSizeArr(setItr) = sum(vocabSizes);
        
        % average mdl score
        vocabSizes(vocabSizes>topXMdl) = topXMdl;
        numberOfSelectedParts = sum(vocabSizes);
        vocabSizes = num2cell(vocabSizes);
        avgMdlScores = cellfun(@(x,y) sum([x(1:y).normMdlScore]), vocabulary, vocabSizes);
        results.mdlScoreArr(setItr) = {avgMdlScores ./ cell2mat(vocabSizes)};
        results.avgMdlScoreArr(setItr) = (sum(avgMdlScores)/numberOfSelectedParts);
        
        %% Time
        load([outputFolder '/trtime.mat']);
        load([outputFolder '/tetime.mat']);
        results.trainTimeArr(setItr) = tr_stop_time;
        results.testTimeArr(setItr) = testInferenceTime;
        
        %% Calculate shareability. 
        % A little tricky, since this depends on the experiment type.
        load([outputFolder '/vb.mat'], 'mainGraph');
        mainGraphSubSets = cellfun(@(x,y) x([x.labelId]<=y), mainGraph, vocabSizes, 'UniformOutput', false);
        levelSizes = cellfun(@(x) numel(x), mainGraphSubSets);
        
        % Assign labels to realizations.
        imageIds = cellfun(@(x) [x.imageId], mainGraphSubSets, 'UniformOutput', false);
        labelIds = cellfun(@(x) [x.labelId], mainGraphSubSets, 'UniformOutput', false);
        load([outputFolder '/vb.mat'], 'trainingFileNames');
        imageLabels = zeros(1, numel(trainingFileNames));
        for fileItr = 1:numel(trainingFileNames)
            [~, fileName, ~] = fileparts(trainingFileNames{fileItr});
            imageLabels(fileItr) = sscanf(fileName, '%d_~%d');
        end
        
        %% Calculate shareability of parts among classes.
        numberOfTotalClasses = numel(unique(imageLabels));
        classLabels = cellfun(@(x) imageLabels(x), imageIds, 'UniformOutput', false);
        shareabilityArr = zeros(numel(levelSizes),1);
        for levelItr = 1:numel(levelSizes)
            levelLabelIds = labelIds{levelItr};
            levelClassLabels = classLabels{levelItr};
            numberOfClasses = 0;
            for partItr = 1:numel(vocabSizes{levelItr})
                numberOfClasses = numberOfClasses + numel(unique(levelClassLabels(levelLabelIds==partItr)));
            end
            shareabilityArr(levelItr) = numberOfClasses / (numberOfTotalClasses^2);
        end
        results.shareabilityArr{setItr} = shareabilityArr;
        results.avgShareabilityArr(setItr) = sum(shareabilityArr .* cell2mat(vocabSizes))/sum(cell2mat(vocabSizes));
        clear vocabulary mainGraph
    end
    save([pwd '/output_ECCV2014/' expType 'results.mat'], 'results');
end

