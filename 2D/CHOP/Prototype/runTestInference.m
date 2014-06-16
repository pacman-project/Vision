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
    global vocabulary, global redundantVocabulary, global modes, global highLevelModes;

    %% ========== Step 1: Run inference for all test images with the learned vocabulary. ==========
    options = SetParameters(datasetName, false);
    
    % override options.
    if ~options.fastInference
%       options.noveltyThr = options.noveltyThr/2;
        options.noveltyThr = 0;
    end
    
    if options.testImages
        %% Step 1.0: Read vocabulary if it exists.
        testFileNames = fuf([options.currentFolder '/input/' datasetName '/test/*' ext], 1, 'detail');
        if exist([options.currentFolder '/output/' datasetName '/vb.mat'], 'file')
            load([options.currentFolder '/output/' datasetName '/vb.mat']);
        else
            display('No vocabulary exists!');
        end
        
        %% Step 1.2: Run inference on each test image.
        totalInferenceTime = 0;
        %% TODO: Parallelize? 
        for testImgItr = 1:size(testFileNames,1) 
            totalInferenceTime = totalInferenceTime + singleTestImage(testFileNames{testImgItr}, options);
        end
        save([options.currentFolder '/output/' datasetName '/tetime.mat'], 'totalInferenceTime');
    end
end