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
            load([options.currentFolder '/output/' datasetName '/vb.mat'], 'modes', 'vocabulary', 'redundantVocabulary', 'categoryNames');
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
        for fileItr = 1:numel(testFileNames)
            categoryLabel = -1; %#ok<NASGU>
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
                chosenCategoryArr = strfind(categoryNames, categoryStr);
                chosenCategoryArr = cellfun(@(x) ~isempty(x), chosenCategoryArr);
                [categoryLabel] = find(chosenCategoryArr, 1, 'first'); %#ok<NASGU>
            end
            save([options.testInferenceFolder '/' fileName '_test.mat'], 'categoryLabel');
        end
        
        %% Step 1.2: Run inference on each test image.
        totalInferenceTime = 0;
        for testImgItr = 1:size(testFileNames,1) 
            totalInferenceTime = totalInferenceTime + singleTestImage(testFileNames{testImgItr}, vocabulary, redundantVocabulary, modes, options);
        end
        save([options.currentFolder '/output/' datasetName '/tetime.mat'], 'totalInferenceTime');
    end
    
    % Close thread pool if opened.
    if options.parallelProcessing
        matlabpool close;
    end
end