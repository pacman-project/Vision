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
   % Read all relevant structures.
   if exist([pwd  '/output/' datasetName '/vb.mat'], 'file')
       load([pwd '/output/' datasetName '/vb.mat'], 'vocabulary', 'allModes', 'modeProbs', 'options', 'categoryNames');
       load([pwd '/output/' datasetName '/distributions.mat'], 'vocabularyDistributions');
   else
       display('No vocabulary exists!');
       totalInferenceTime = 0;
       return;
   end
   options.isTraining = false;
   
       %% ========== PATH FOLDER ADDITION ==========
    w = warning('off', 'all');
    addpath(genpath([options.currentFolder '/utilities']));
    addpath(genpath([options.currentFolder '/demo']));
    addpath(genpath([options.currentFolder '/graphTools']));
    addpath(genpath([options.currentFolder '/vocabLearning']));
    addpath(genpath([options.currentFolder '/inference']));
    addpath(genpath([options.currentFolder '/categorization']));
    warning(w);
   
    % Remove output folder.
    if exist(options.testInferenceFolder, 'dir')
         rmdir(options.testInferenceFolder, 's');
    end
    mkdir(options.testInferenceFolder);
    
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
        catch %#ok<CTCH>
           display('No files under /input/test folder. Please put your test images under this folder.'); 
           totalInferenceTime = 0;
           return;
        end
    
        % Create the folder structure required.
        createFolders(options);
        
        % Learn the category of the images, if they are structured in
        % manner that allows measuring performance.
        categoryArrIdx = zeros(numel(testFileNames),1);
        for fileItr = 1:numel(testFileNames)
            categoryLabel = 1;
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
        
        % For some weird reason, Matlab workers cannot access variables
        % read from the file. They have to be used in the code. Here's my
        % workaround: 
        vocabulary = vocabulary; %#ok<ASGSL,NODEF>
        allModes = allModes; %#ok<ASGSL,NODEF>
        modeProbs = modeProbs; %#ok<ASGSL,NODEF>
        categoryNames = categoryNames; %#ok<ASGSL>
        
        %% Step 1.2: Run inference on each test image.
        startTime = tic;
        for testImgItr = 1:size(testFileNames,1) 
            [~, testFileName, ~] = fileparts(testFileNames{testImgItr});
            display(['Processing ' testFileName '...']);
            singleTestImage(testFileNames{testImgItr}, vocabulary, vocabularyDistributions, categoryNames{categoryArrIdx(testImgItr)}, options); 
        end
        totalInferenceTime = toc(startTime);
        save([options.currentFolder '/output/' datasetName '/tetime.mat'], 'totalInferenceTime');
    end
    
    % Close thread pool if opened.
    if options.parallelProcessing
        matlabpool close;
    end
end