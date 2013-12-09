%> Name: visualizeImages
%>
%> Description: Given the current graph and vocabulary, this function
%> intends to visualize each image by putting the feature patches in
%> estimated positions.
%>
%> @param fileList Image list to work on.
%> @param vocabCompCount The number of words in current vocabulary
%> level.
%> @param graphLevel The graph level which includes specific graph of the
%> each image. Assumed to be ordered by image id.
%> @param levelItr The level index.
%> @param options Program options.
%> @param datasetname The name of the dataset.
%> @param type 'test' if the image list consist of test images, 'train' if 
%> training.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 09.12.2013
function [ ] = visualizeImages( fileList, vocabCompCount, graphLevel, levelItr, options, datasetName, type )
    outputDir = [options.outputFolder '/reconstruction/' type];
    if ~exist(outputDir, 'dir')
       mkdir(outputDir); 
    end
    processedDir = options.processedFolder;
    
    %% Read vocabulary into separate masks to avoid re-reading of images.
    vocabMasks = cell(vocabCompCount,1);
    vocabMaskDir = [options.currentFolder '/debug/' datasetName '/level' num2str(levelItr) '/reconstruction'];
    for fileItr = 1:vocabCompCount
        vocabMasks{fileItr} = imread([vocabMaskDir '/' num2str(fileItr) '.png']);
    end
    
    for fileItr = 1:numel(fileList)
        %% Learn the size of the original image, and allocate space for new mask.
        [~, fileName, ~] = fileparts(fileList{fileItr});
        img = imread([processedDir '/' fileName '.png']);
        reconstructedMask = zeros(size(img,1), size(img,2), 'uint8');
        
        %% Go over each node and write its mask to the image's reconstruction.
        nodes = graphLevel(:,[graphLevel.imageId] == fileItr);
        for nodeItr = 1:numel(nodes)
            % Read the mask here, and crop it if necessary.
            nodeMask = vocabMasks{nodes(nodeItr).labelId};
            
            halfSize = (size(nodeMask)-1)/2;
            minX = halfSize(1);
            maxX = halfSize(1);
            minY = halfSize(2);
            maxY = halfSize(2);
            imageSize = size(reconstructedMask);
            
            if imageSize(1)<nodes(nodeItr).position(1)+maxX;
                maxX = imageSize(1) - nodes(nodeItr).position(1);
            end
            if 1>(nodes(nodeItr).position(1)-minX);
                minX = nodes(nodeItr).position(1) - 1;
            end
            if imageSize(2)<(nodes(nodeItr).position(2)+maxY);
                maxY = imageSize(2) - nodes(nodeItr).position(2);
            end
            if 1>(nodes(nodeItr).position(2)-minY);
                minY = nodes(nodeItr).position(2) - 1;
            end
            
%             reconstructedMask((nodes(nodeItr).position(1)-minX):(nodes(nodeItr).position(1)+maxX), ...
%                 (nodes(nodeItr).position(2)-minY):(nodes(nodeItr).position(2)+maxY)) = ...
%                 max(reconstructedMask((nodes(nodeItr).position(1)-minX):(nodes(nodeItr).position(1)+maxX), ...
%                 (nodes(nodeItr).position(2)-minY):(nodes(nodeItr).position(2)+maxY)), ...
%                 nodeMask(((halfSize(1)-minX)+1):(end-(halfSize(1)-maxX)), ...
%                 ((halfSize(2)-minY)+1):(end-(halfSize(2)-maxY))));
            reconstructedMask((nodes(nodeItr).position(1)-minX):(nodes(nodeItr).position(1)+maxX), ...
                (nodes(nodeItr).position(2)-minY):(nodes(nodeItr).position(2)+maxY)) = ...
                nodeMask(((halfSize(1)-minX)+1):(end-(halfSize(1)-maxX)), ...
                ((halfSize(2)-minY)+1):(end-(halfSize(2)-maxY)));
        end
        
        %% Write the reconstructed mask to the output.
        imwrite(reconstructedMask, [outputDir, '/' fileName '_level' num2str(levelItr) '.png']);
    end
end

