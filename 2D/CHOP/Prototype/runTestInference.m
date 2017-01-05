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
       %% ========== PATH FOLDER ADDITION ==========
    w = warning('off', 'all');
    addpath(genpath([pwd '/utilities']));
    addpath(genpath([pwd '/demo']));
    addpath(genpath([pwd '/graphTools']));
    addpath(genpath([pwd '/vocabLearning']));
    addpath(genpath([pwd '/inference']));
    addpath(genpath([pwd '/categorization']));
    warning(w);
    
   % Read all relevant structures.
   if exist([pwd  '/output/' datasetName '/vb.mat'], 'file')
       load([pwd '/output/' datasetName '/vb.mat'], 'vocabulary', 'options', 'categoryNames');
       load([pwd '/output/' datasetName '/distributions.mat'], 'vocabularyDistributions');
       load([pwd '/output/' datasetName '/export.mat'], 'exportArr', 'trainingFileNames', 'activationArr');
   else
       display('No vocabulary exists!');
       totalInferenceTime = 0;
       return;
   end
   options.isTraining = false;
   options.testDebug = false;
   
    % Remove output folder.
    if exist(options.testInferenceFolder, 'dir')
         rmdir(options.testInferenceFolder, 's');
    end
    mkdir(options.testInferenceFolder);
% %     
    % Open threads for parallel processing.
    if options.parallelProcessing && usejava('jvm')
        if exist('parpool', 'file')
            p = gcp('nocreate');
            if ~isempty(p)
                delete(p);
            end
            parpool(options.numberOfThreads);
        else
            s = matlabpool('size');
            if s>0
               matlabpool close; 
            end
            matlabpool('open', options.numberOfThreads);
        end
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
            strLength = numel([options.currentFolder '/input/' datasetName '/test/']);
            fileNameLength = numel(ext) + numel(fileName) + 1; % 1 for '/' (folder delimiter)
            if numel(fullName) >= strLength + fileNameLength
                categoryStr = fullName(1, (strLength+1):(end-fileNameLength));
                categoryStrSepIdx = union(strfind(categoryStr, '/'), strfind(categoryStr, '\'));
                if ~isempty(categoryStrSepIdx)
                    categoryStr = categoryStr(:, 1:(categoryStrSepIdx(1)-1));
                end
                chosenCategoryArr = cellfun(@(x) strcmp(x, categoryStr), categoryNames);
                categoryLabel = find(chosenCategoryArr, 1, 'first');
            end
            if ~isempty(categoryLabel)
                categoryArrIdx(fileItr) = categoryLabel;
            end
        end
        
        % Get separation. 
        if strfind(trainingFileNames{1}, '/')
             sep = '/';
        else
             sep = '\';
        end
        
        %% For fast processing, we put vocabulary children in arrays.
        vocabularyChildren = cell(numel(vocabulary), 1);
        uniqueVocabularyChildren = cell(numel(vocabulary), 1);
        for levelItr = 1:numel(vocabulary)
             if isempty(vocabulary{levelItr})
                  break;
             end
             vocabLevelChildren = {vocabulary{levelItr}.children};
             maxCount = max(cellfun(@(x) numel(x), vocabLevelChildren));
             vocabLevelChildren = cellfun(@(x) [x, zeros(1, maxCount - numel(x))], vocabLevelChildren, 'UniformOutput', false);
             vocabLevelChildren = cat(1, vocabLevelChildren{:});
             vocabularyChildren{levelItr} = vocabLevelChildren;
             vocabLevelChildren = unique(vocabLevelChildren, 'rows');
             uniqueVocabularyChildren{levelItr} = vocabLevelChildren;
        end
        
        % For some weird reason, Matlab workers cannot access variables
        % read from the file. They have to be used in the code. Here's my
        % workaround: 
        vocabulary = vocabulary; %#ok<ASGSL,NODEF>
        vocabularyDistributions = vocabularyDistributions; %#ok<ASGSL,NODEF>
        categoryNames = categoryNames; %#ok<ASGSL>
        
        %% Step 1.2: Run inference on each test image.
        startTime = tic;
        confMatrix = zeros(numel(categoryNames),numel(categoryNames));
        categoryInfo = cell(size(testFileNames,1) ,1);
        allExportArr = exportArr;
        allActivationArr = activationArr;
        top1CorrectArr = zeros(size(testFileNames,1),1) > 0;
        top3CorrectArr = zeros(size(testFileNames,1),1) > 0;
        top5CorrectArr = zeros(size(testFileNames,1),1) > 0;
        parfor testImgItr = 1:size(testFileNames,1) 
    %    for testImgItr = 1:5
            [~, testFileName, ~] = fileparts(testFileNames{testImgItr});
            display(['Processing ' testFileName '...']);
            [categoryLabel, predictedCategory, top3Correct, top5Correct] = singleTestImage(testFileNames{testImgItr}, vocabulary, vocabularyDistributions, uniqueVocabularyChildren, vocabularyChildren, categoryArrIdx(testImgItr), categoryNames{categoryArrIdx(testImgItr)}, options);
           
%             imageId = find(cellfun(@(x) ~isempty(x), strfind(trainingFileNames, [sep testFileName ext])));
%             if ~isempty(imageId)
%                  imageExportArr = allExportArr(allExportArr(:,5) == imageId,:);
%                  imageActivationArr = allActivationArr(allExportArr(:,5) == imageId, :);
%             end
    %        load([options.testInferenceFolder '/' categoryNames{categoryArrIdx(testImgItr)} '_' testFileName '_test.mat'], 'exportArr', 'activationArr');
    %        testExportArr = exportArr;
%            testActivationArr = activationArr;
%            if ~isempty(imageExportArr)
%                 [sortedImageExportArr, ~] = sortrows(imageExportArr);
%                 sortedImageExportArr = sortedImageExportArr(:, 1:4);
%                 [sortedTestExportArr, ~] = sortrows(testExportArr);
%                 sortedTestExportArr = sortedTestExportArr(:, 1:4);
%                 if ~isequal(sortedImageExportArr, sortedTestExportArr)
%                      display(['Image ' num2str(testImgItr) ' failed.']);
%  %               elseif ~isequal(imageActivationArr(idx1), testActivationArr(idx2))
%    %                  display(['Image ' num2str(testImgItr) ' failed.']);
%                 end
%            end
%            
    %        maxLevel = max(testExportArr(:,4));
            categoryInfo(testImgItr) = {[categoryLabel, predictedCategory]};
            top1CorrectArr(testImgItr) = categoryLabel == predictedCategory;
            top3CorrectArr(testImgItr) = top3Correct;
            top5CorrectArr(testImgItr) = top5Correct;
        end
        totalInferenceTime = toc(startTime);
        categoryInfo = cat(1, categoryInfo{:});
        for itr = 1:size(categoryInfo,1)
            confMatrix(categoryInfo(itr,1), categoryInfo(itr,2)) = confMatrix(categoryInfo(itr,1), categoryInfo(itr,2)) + 1;
        end
        perf = sum(confMatrix(1:(size(confMatrix,1)+1):(size(confMatrix,1))^2)) / size(testFileNames,1); %#ok<NASGU>
        top1Accuracy = nnz(top1CorrectArr) / numel(top1CorrectArr);
        top3Accuracy = nnz(top3CorrectArr) / numel(top3CorrectArr);
        top5Accuracy = nnz(top5CorrectArr) / numel(top5CorrectArr);
        save([options.currentFolder '/output/' datasetName '/classification.mat'], 'confMatrix', 'perf', 'top1CorrectArr', 'top1Accuracy', 'top3CorrectArr', 'top3Accuracy', 'top5CorrectArr', 'top5Accuracy');
        save([options.currentFolder '/output/' datasetName '/tetime.mat'], 'totalInferenceTime');
    end
    
    % Close thread pool if opened.
    if options.parallelProcessing
        if ~exist('parpool', 'file')
            matlabpool close;
        end
    end
end

function [predictedCategory] = parload(pathFile, varName)
    predictedCategory = -1;
    load(pathFile,varName);
end