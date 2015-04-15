%> Name: runTestInference
%>
%> Description: Assuming the unsupervised vocabulary has been learned up to
%> a certain level, this function runs inference on test images to find
%> realizations of parts in vocabulary. This process is unsupervised, and
%> comes to an end in the last level in the vocabulary.
%>
%> @param datasetName Name of the dataset to work on. 
%> @param ext file extension.
%> 
%> @retval options Program options.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 01.02.2014
%> Ver 1.1 on 12.01.2014 Timing is added by Mete
function [ totalInferenceTime ] = runTestInference( datasetName, ext )
    %% ========== Step 1: Run inference for all test images with the learned vocabulary. ==========
    options = SetParameters(datasetName, false);
    options.currentFolder = options.currentFolder;
    
    %% Learn edge-based distance matrix once and for all.
    [edgeIdMatrix, edgeDistanceMatrix, edgeCoords] = findEdgeDistanceMatrix(options.edgeQuantize);
    options.edgeIdMatrix = edgeIdMatrix;
    options.edgeDistanceMatrix = edgeDistanceMatrix;
    options.edgeCoords = edgeCoords;
    clear edgeIdMatrix edgeDistanceMatrix edgeCoords;
    
    % Open threads for parallel processing.
    if options.parallelProcessing
        s = matlabpool('size');
        if s>0
           matlabpool close; 
        end
        matlabpool('open', options.numberOfThreads);
    end
    
    if options.testImages
        %% Step 1.0: Read vocabulary if it exists.
        try
            testFileNames = fuf([options.currentFolder '/input/' datasetName '/test/*' ext], 1, 'detail');
        catch
           display('No files under /input/test folder. Please put your test images under this folder.'); 
           totalInferenceTime = 0;
           return;
        end
    
        if exist([options.currentFolder '/output/' datasetName '/vb.mat'], 'file')
            load([options.currentFolder '/output/' datasetName '/vb.mat'], 'vocabulary', 'distanceMatrices', 'graphLevelIndices', 'categoryNames', 'optimalThresholds');
            load([options.currentFolder '/output/' datasetName '/export.mat']);
        else
            display('No vocabulary exists!');
            totalInferenceTime = 0;
            return;
        end
        
        % Create the folder structure required.
        createFolders(options);
        
        % Learn the category of the images, if they are structured in
        % manner that allows measuring performance.
        categoryArrIdx = zeros(numel(testFileNames),1);
        for fileItr = 1:numel(testFileNames)
            categoryLabel = -1;
            fullName = testFileNames{fileItr};
            [~, fileName, ext] = fileparts(testFileNames{fileItr});
            strLength = numel([options.currentFolder '/input/' options.datasetName '/test/']);
            fileNameLength = numel(ext) + numel(fileName) + 1; % 1 for '/' (folder delimiter)
            if numel(fullName) >= strLength + fileNameLength
                categoryStr = fullName(1, (strLength+1):(end-fileNameLength));
                categoryStrSepIdx = strfind(categoryStr, '/');
                if ~isempty(categoryStrSepIdx)
                    categoryStr = categoryStr(:, 1:(categoryStrSepIdx(1)-1));
                end
                chosenCategoryArr = cellfun(@(x) strcmp(x, categoryStr), categoryNames);
                categoryLabel = find(chosenCategoryArr, 1, 'first');
            end
            if ~isempty(categoryLabel)
                categoryArrIdx(fileItr) = categoryLabel;
            end
            save([options.testInferenceFolder '/' categoryNames{categoryLabel} '_' fileName '_test.mat'], 'categoryLabel');
        end
        
        %% We have modified the ranking in vocabulary based on frequency, after learning them using MDL. 
        % Now, we revert back to the original ranking.
        if numel(vocabulary)>1 %#ok<NODEF>
            for levelItr = 2:numel(vocabulary)
                vocabLevel = vocabulary{levelItr};
                arrToSort = [vocabLevel.orgRank];
                [~, sortIdx] = sort(arrToSort, 'ascend');
                vocabLevel = vocabLevel(sortIdx);
                % Assign labels.
                vocabulary(levelItr) = {vocabLevel}; %#ok<AGROW>
            end
        end
        
        % For some weird reason, Matlab workers cannot access variables
        % read from the file. They have to be used in the code. Here's my
        % workaround: 
        distanceMatrices = distanceMatrices; %#ok<NODEF,ASGSL>
        optimalThresholds = optimalThresholds; %#ok<NODEF,ASGSL>
        categoryNames = categoryNames; %#ok<ASGSL>
        
        %% Step 1.2: Run inference on each test image.
        startTime = tic;
        parfor testImgItr = 1:size(testFileNames,1) 
            [~, testFileName, ~] = fileparts(testFileNames{testImgItr});
            display(['Processing ' testFileName '...']);
            singleTestImage(testFileNames{testImgItr}, vocabulary, distanceMatrices, categoryNames{categoryArrIdx(testImgItr)}, optimalThresholds, options); 
        end
        totalInferenceTime = toc(startTime);
        save([options.currentFolder '/output/' datasetName '/tetime.mat'], 'totalInferenceTime');
    end
    
    % Close thread pool if opened.
    if options.parallelProcessing
        matlabpool close;
    end
end