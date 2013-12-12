%> Name: visualizeImages
%>
%> Description: Given the current graph and vocabulary, this function
%> intends to visualize each image by putting the feature patches in
%> estimated positions.
%>
%> @param fileList Image list to work on.
%> @param mainGraph The main graph which includes specific graph of the
%> each image at each level. Each level is assumed to be ordered by image id.
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
function [ ] = visualizeImages( fileList, mainGraph, levelItr, options, datasetName, type )
    outputDir = [options.outputFolder '/reconstruction/' type];
    if ~exist(outputDir, 'dir')
       mkdir(outputDir); 
    end
    processedDir = options.processedFolder;
    
    firstLevel = mainGraph{1};
    graphLevel = mainGraph{levelItr};
    
    %% Read level 1 into separate masks.
    vocabMasks = cell(options.numberOfFilters,1);
    vocabMaskDir = [options.currentFolder '/debug/' datasetName '/level1/reconstruction'];
    for fileItr = 1:options.numberOfFilters
        vocabMasks{fileItr} = imread([vocabMaskDir '/' num2str(fileItr) '.png']);
    end
    
    for fileItr = 1:numel(fileList)
        %% Learn the size of the original image, and allocate space for new mask.
        [~, fileName, ~] = fileparts(fileList{fileItr});
        img = imread([processedDir '/' fileName '.png']);
        reconstructedMask = zeros(size(img,1), size(img,2), 'uint8');
        labeledReconstructedMask = zeros(size(img,1), size(img,2));
        
        %% Go over each leaf in the graph and write its mask to the reconstructed image.
        nodes = graphLevel(:,[graphLevel.imageId] == fileItr);
        
        for nodeItr = 1:numel(nodes)
            leafNodes = nodes(nodeItr).leafNodes;
            for leafNodeItr = 1:numel(leafNodes)
                % Read the mask here, and crop it if necessary.
                nodeMask = vocabMasks{firstLevel(leafNodes(leafNodeItr)).labelId};
                position = firstLevel(leafNodes(leafNodeItr)).position;
                
                % If the part has already been added, move on.
                if reconstructedMask(position(1), position(2))>0
                   continue; 
                end

                halfSize = (size(nodeMask)-1)/2;
                minX = halfSize(1);
                maxX = halfSize(1);
                minY = halfSize(2);
                maxY = halfSize(2);
                imageSize = size(reconstructedMask);

                if imageSize(1)<position(1)+maxX;
                    maxX = imageSize(1) - position(1);
                end
                if 1>(position(1)-minX);
                    minX = position(1) - 1;
                end
                if imageSize(2)<(position(2)+maxY);
                    maxY = imageSize(2) - position(2);
                end
                if 1>(position(2)-minY);
                    minY = position(2) - 1;
                end
                 reconstructedMask((position(1)-minX):(position(1)+maxX), ...
                     (position(2)-minY):(position(2)+maxY)) = ... 
                     max(reconstructedMask((position(1)-minX):(position(1)+maxX), ...
                     (position(2)-minY):(position(2)+maxY)), ...
                     nodeMask(((halfSize(1)-minX)+1):(end-(halfSize(1)-maxX)), ...
                     ((halfSize(2)-minY)+1):(end-(halfSize(2)-maxY))));
                 
                 % First, write the node label to the labeled mask.
                 reconstructedPatch = labeledReconstructedMask((position(1)-minX):(position(1)+maxX), ...
                     (position(2)-minY):(position(2)+maxY));
                 reconstructedPatch(nodeMask(((halfSize(1)-minX)+1):(end-(halfSize(1)-maxX)), ...
                     ((halfSize(2)-minY)+1):(end-(halfSize(2)-maxY))) > 10) = nodeItr;
                 labeledReconstructedMask((position(1)-minX):(position(1)+maxX), ...
                     (position(2)-minY):(position(2)+maxY)) = reconstructedPatch;
                 
                 % Write the response to the response mask.
                 reconstructedMask((position(1)-minX):(position(1)+maxX), ...
                     (position(2)-minY):(position(2)+maxY)) = ...
                     max(reconstructedMask((position(1)-minX):(position(1)+maxX), ...
                     (position(2)-minY):(position(2)+maxY)), ...
                     nodeMask(((halfSize(1)-minX)+1):(end-(halfSize(1)-maxX)), ...
                     ((halfSize(2)-minY)+1):(end-(halfSize(2)-maxY))));
    %            reconstructedMask((position(1)-minX):(position(1)+maxX), ...
    %                (position(2)-minY):(position(2)+maxY)) = ...
    %                nodeMask(((halfSize(1)-minX)+1):(end-(halfSize(1)-maxX)), ...
    %                ((halfSize(2)-minY)+1):(end-(halfSize(2)-maxY)));
            end
        end
        
        %% Write the reconstructed mask to the output.
        % Add some random colors to make each composition look different.
        rgbImg = label2rgb(labeledReconstructedMask, 'jet', 'k', 'shuffle');
        rgbImg(:,:,1) = uint8(double(rgbImg(:,:,1)) .* (double(reconstructedMask)/255));
        rgbImg(:,:,2) = uint8(double(rgbImg(:,:,2)) .* (double(reconstructedMask)/255));
        rgbImg(:,:,3) = uint8(double(rgbImg(:,:,3)) .* (double(reconstructedMask)/255));
        
        if levelItr>1
            % Read first level to fill the voids left by missing filters.
            firstLevelImg = imread([outputDir, '/' fileName '_level1.png']);
            firstLevelMask=max(firstLevelImg, [], 3);
   %         firstLevelMask = uint8(double(firstLevelMask)/double(max(max(firstLevelMask))));
            firstLevelMask = (firstLevelMask .* uint8(labeledReconstructedMask==0));
            firstLevelImg(:,:,1) = firstLevelMask;
            firstLevelImg(:,:,2) = firstLevelMask;
            firstLevelImg(:,:,3) = firstLevelMask;

            % Combine both and write to output.
            imwrite(rgbImg + firstLevelImg, [outputDir, '/' fileName '_level' num2str(levelItr) '.png']);
        else
            % Write to output.
            imwrite(rgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) '.png']);
        end
    end
end

