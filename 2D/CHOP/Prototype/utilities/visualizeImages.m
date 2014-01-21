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
function [ ] = visualizeImages( fileList, mainGraph, levelItr, options, ~, type )
    outputDir = [options.outputFolder '/reconstruction/' type];
    originalDir = [options.outputFolder '/original/'];
    if ~exist(outputDir, 'dir')
       mkdir(outputDir); 
    end
    processedDir = options.processedFolder;
    
    firstLevel = mainGraph{1};
    graphLevel = mainGraph{levelItr};
    
    %% Read level 1 into separate masks.
    numberOfFilters = getNumberOfFilters(options);
    vocabMasks = cell(numberOfFilters,1);
    vocabMaskDir = [options.currentFolder '/filters/' options.filterType];
    for fileItr = 1:numberOfFilters
        vocabMasks{fileItr} = imread([vocabMaskDir '/filt' num2str(fileItr) '.png']);
    end
    
    for fileItr = 1:numel(fileList)
        %% Learn the size of the original image, and allocate space for new mask.
        [~, fileName, ~] = fileparts(fileList{fileItr});
        img = imread([processedDir '/' fileName '.png']);
        originalImg = imread([originalDir fileName '.png']);
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
                imageSize = size(reconstructedMask);

                %% If patch is out of bounds, do nothing.
                if imageSize(1) < position(1)+halfSize(1) || ...
                        1 > position(1)-halfSize(1) || ...
                        imageSize(2) < (position(2)+halfSize(2)) ||...
                        1 > (position(2)-halfSize(2))
                    continue;
                end
                
                reconstructedMask((position(1)-halfSize(1)):(position(1)+halfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+halfSize(2))) = ... 
                     max(reconstructedMask((position(1)-halfSize(1)):(position(1)+halfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+halfSize(2))), ...
                     nodeMask);
                 
                 % First, write the node label to the labeled mask.
                 reconstructedPatch = labeledReconstructedMask((position(1)-halfSize(1)):(position(1)+halfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+halfSize(2)));
                 reconstructedPatch(nodeMask > 10) = nodeItr;
                 labeledReconstructedMask((position(1)-halfSize(1)):(position(1)+halfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+halfSize(2))) = reconstructedPatch;
            end
        end
        
        %% Write the reconstructed mask to the output.
        % Add some random colors to make each composition look different, 
        % and overlay the gabors with the original image.
        rgbImg = label2rgb(labeledReconstructedMask, 'jet', 'k', 'shuffle');
        %% For level 1, write the original image too.
        if size(originalImg,3) == 1
           oldOriginalImg = originalImg;
           originalImg = zeros(size(originalImg,1), size(originalImg,2), 3, 'uint8');
           for bandItr = 1:3
              originalImg(:,:,bandItr) = oldOriginalImg; 
           end
        end
        
        for bandItr = 1:size(rgbImg,3)
            rgbImg(:,:,bandItr) = uint8(double(rgbImg(:,:,bandItr)) .* (double(reconstructedMask)/255)) + ...
            originalImg(:,:,bandItr) .* uint8(reconstructedMask==0);
        end 
       
        if levelItr>1
            % Read first level to fill the voids left by missing filters.
            firstLevelMask = imread([outputDir, '/' fileName '_level1clean.png']);
            firstLevelImg = zeros(size(firstLevelMask,1), size(firstLevelMask,2), size(rgbImg,3), 'uint8');
            for bandItr = 1:size(rgbImg,3)
                firstLevelImg(:,:,bandItr) = firstLevelMask;
            end

            % Combine both and write to output.
            imwrite(rgbImg + firstLevelImg, [outputDir, '/' fileName '_level' num2str(levelItr) '.png']);
        else
            imwrite(rgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) '.png']);
            imwrite(reconstructedMask, [outputDir, '/' fileName '_level' num2str(levelItr) 'clean.png']);
        end
    end
end