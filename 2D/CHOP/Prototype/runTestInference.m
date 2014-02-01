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
%> Ver 1.0 on 10.01.2014
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
        preparePreDefinedFiles(options.preDefinedFolder, vocabulary);
        
        %% Step 1.2: Run inference on each test image.
        for testImgItr = 1:size(testFileNames,1)
            singleTestImage(testFileNames{testImgItr}, options, options.currentFolder);
        end
    end
end

