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
       load([pwd '/output/' datasetName '/vb.mat'], 'vocabulary', 'allModes', 'modeProbs', 'options', 'categoryNames');
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
            if ~exist([options.testInferenceFolder '/' categoryNames{categoryLabel} '_' fileName '_test.mat'], 'file')
                 save([options.testInferenceFolder '/' categoryNames{categoryLabel} '_' fileName '_test.mat'], 'categoryLabel');
            else
                 save([options.testInferenceFolder '/' categoryNames{categoryLabel} '_' fileName '_test.mat'], 'categoryLabel', '-append');
            end
        end
        
        % Get separation. 
        if strfind(trainingFileNames{1}, '/')
             sep = '/';
        else
             sep = '\';
        end
        
        %% Finally, we modify edge info so we can match nodes faster.
        for vocabLevelItr = 2:numel(vocabulary)
             vocabLevel = vocabulary{vocabLevelItr};
             
             if isempty(vocabLevel)
                  break;
             end
             
             modes = allModes{vocabLevelItr-1};
             [uniqueModes, IC, ~] = unique(modes(:,1:2), 'rows', 'R2012a');
             IC = cat(1, IC, size(modes,1)+1);
             for nodeItr = 1:numel(vocabLevel)
                  vocabNode = vocabLevel(nodeItr);
                  curAdjInfo = vocabNode.adjInfo;
                  curChildren = vocabNode.children;
                  for edgeItr = 1:size(curAdjInfo,1)
                        relevantModeId = find(uniqueModes(:,1) == curChildren(1) & uniqueModes(:,2) == curChildren(1+edgeItr));
                        subModeIds = modes(IC(relevantModeId):(IC(relevantModeId+1)-1), 3);
                        subModeId = find(subModeIds == curAdjInfo(edgeItr,3));
                        if ~isempty(subModeId)
                             curAdjInfo(edgeItr,3) = subModeId;
                        else
                             display('Error in edge assignment!');
                        end
                  end
                  vocabNode.adjInfo = curAdjInfo;
                  vocabLevel(nodeItr) = vocabNode;
             end
             vocabulary{vocabLevelItr} = vocabLevel;
        end
        
        %% Create mode indices for fast mode indexing.
        
        % For some weird reason, Matlab workers cannot access variables
        % read from the file. They have to be used in the code. Here's my
        % workaround: 
        vocabulary = vocabulary; %#ok<ASGSL,NODEF>
        vocabularyDistributions = vocabularyDistributions; %#ok<ASGSL,NODEF>
        allModes = allModes; %#ok<ASGSL,NODEF>
        modeProbs = modeProbs; %#ok<ASGSL,NODEF>
        categoryNames = categoryNames; %#ok<ASGSL>
        
        %% Step 1.2: Run inference on each test image.
        startTime = tic;
        confMatrix = zeros(numel(categoryNames),numel(categoryNames));
        categoryInfo = cell(size(testFileNames,1) ,1);
        parfor testImgItr = 1:size(testFileNames,1) 
    %    for testImgItr = 1:5
            [~, testFileName, ~] = fileparts(testFileNames{testImgItr});
            display(['Processing ' testFileName '...']);
            [categoryLabel, predictedCategory] = singleTestImage(testFileNames{testImgItr}, vocabulary, vocabularyDistributions, allModes, modeProbs, categoryNames{categoryArrIdx(testImgItr)}, options);
           
%             imageId = find(cellfun(@(x) ~isempty(x), strfind(trainingFileNames, [sep testFileName ext])));
%             if ~isempty(imageId)
%                  imageExportArr = exportArr(exportArr(:,5) == imageId,:);
%                  imageActivationArr = activationArr(exportArr(:,5) == imageId, :);
%             end
%            load([options.testInferenceFolder '/' categoryNames{categoryArrIdx(testImgItr)} '_' testFileName '_test.mat'], 'exportArr', 'activationArr');
%            testExportArr = exportArr;
%            testActivationArr = activationArr;
%            if ~isempty(imageExportArr)
%                 [sortedImageExportArr, idx1] = sortrows(imageExportArr);
%                 idx1 = idx1(sortedImageExportArr(:,4) > 1);
%                 [sortedTestExportArr, idx2] = sortrows(testExportArr);
%                 idx2 = idx2(sortedTestExportArr(:,4) > 1);
%                 if ~isequal(sortedImageExportArr(:,1:4), sortedTestExportArr(:,1:4))
%                      display(['Image ' num2str(testImgItr) ' failed.']);
%                 elseif ~isequal(imageActivationArr(idx1), testActivationArr(idx2))
%                      display(['Image ' num2str(testImgItr) ' failed.']);
%                 end
%            end
           
            categoryInfo(testImgItr) = {[categoryLabel, predictedCategory]};
        end
        totalInferenceTime = toc(startTime);
        categoryInfo = cat(1, categoryInfo{:});
        for itr = 1:size(categoryInfo,1)
            confMatrix(categoryInfo(itr,1), categoryInfo(itr,2)) = confMatrix(categoryInfo(itr,1), categoryInfo(itr,2)) + 1;
        end
        perf = sum(confMatrix(1:(size(confMatrix,1)+1):(size(confMatrix,1))^2)) / size(testFileNames,1); %#ok<NASGU>
        save([options.currentFolder '/output/' datasetName '/classification.mat'], 'confMatrix', 'perf');
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