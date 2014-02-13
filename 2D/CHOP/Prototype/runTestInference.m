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
function [ ] = runTestInference( datasetName, ext )
    %% ========== Step 1: Run inference for all test images with the learned vocabulary. ==========
    options = SetParameters(datasetName);
    if options.testImages
        %% Step 1.0: Read vocabulary if it exists.
        testFileNames = fuf([pwd '/input/' datasetName '/vocab/*' ext], 1, 'detail');
        if exist([options.currentFolder '/output/' datasetName '/' datasetName '_vb.mat'], 'file')
            load([options.currentFolder '/output/' datasetName '/' datasetName '_vb.mat']);
        else
            display('No vocabulary exists!');
        end
        
        %% Step 1.1: Create files for pre-defined substructures ( compositions from voc. at each level)
        if strcmp(options.subdue.implementation, 'exe');
            preparePreDefinedFiles(options.preDefinedFolder, vocabulary);
        end
        
        %% Step 1.2: Run inference on each test image.
        for testImgItr = 1:size(testFileNames,1)
         te_s_time=tic;   
    
            singleTestImage(testFileNames{testImgItr}, options, options.currentFolder);
         test_stop_time(testImgItr)=toc(te_s_time);
        end
        save([options.currentFolder '/output/' datasetName '/' datasetName '_tetime.mat'], 'test_stop_time');
    end
end