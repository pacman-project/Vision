%> Name: visualizeImages
%>
%> Description: Given the current graph and vocabulary, this function
%> intends to visualize each image by putting the feature patches in
%> estimated positions.
%>
%> @param fileList Image list to work on.
%> @param vocabLevel Current vocabulary level.
%> @param graphLevel The object graph which includes specific graph of the
%> each image at each level. Each level is assumed to be ordered by image id.
%> @param levelItr The level index.
%> @param options Program options.
%> @param type 'test' if the image list consist of test images, 'train' if 
%> training.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 09.12.2013
function [ ] = visualizeImages( fileList, ~, graphLevel, leafNodes, levelItr, options, type )
    outputTempDir = [options.outputFolder '/reconstruction/' type];
    backgroundDir = [options.outputFolder '/reconstruction/' type '/' options.backgroundClass];
    processedFolder = options.processedFolder;   
    reconstructionType = options.reconstructionType;
    
    filter1 = options.filters{1};
    filtBandCount = size(filter1,3);
    
    %% Depending on the reconstruction type, we read masks to put on correct positions in images.
    if strcmp(reconstructionType, 'true')
        % Read masks of compositions in current level.
        usedLevel = levelItr;
    else
        % Read level 1 into separate masks.
        usedLevel = 1;
    end
    
    %% Read vocabulary masks.
    vocabMasks = cell(numel(options.filters),1);
    if usedLevel == 1
        for fileItr = 1:numel(options.filters)
            filter1 = options.filters{fileItr};
            filter1 = uint8(255*(filter1 - min(min(min(filter1))))/(max(max(max(filter1))) - min(min(min(filter1)))));
            vocabMasks{fileItr} = filter1;
        end
    else
        vocabMaskList = fuf([options.currentFolder '/debug/' options.datasetName '/level' num2str(usedLevel) '/reconstruction/*.png'], 1, 'detail');
        minusList = fuf([options.currentFolder '/debug/' options.datasetName '/level' num2str(usedLevel) '/reconstruction/*_comp.png'], 1, 'detail');
        numberOfMasks = numel(vocabMaskList) - numel(minusList);
        for fileItr = 1:numberOfMasks
            vocabMasks{fileItr} = imread([options.currentFolder '/debug/' options.datasetName ...
                '/level' num2str(usedLevel) '/reconstruction/' num2str(fileItr) '.png']);
        end
    end
    
    %% Put realizations into distinct sets so that each image has its own nodes in a set.
    imageIds = [graphLevel.imageId]';
    imageNodeSets = cell(numel(fileList),1);
    for fileItr = 1:numel(fileList)
        imageNodeSets(fileItr) = {graphLevel(imageIds == fileItr)};
    end
    
    %% Go over the list of images and run reconstruction.
    parfor fileItr = 1:numel(fileList)
        nodeOffset = numel(find(imageIds<fileItr));
        %% Learn the size of the original image, and allocate space for new mask.
        [~, fileName, ~] = fileparts(fileList{fileItr});
        img = imread([processedFolder '/' fileName '.png']);
        actualImg = img;
        
        % If original image's band count is less than 3, duplicate
        % bands to have a 3-band image.
        if size(actualImg,3) == 1
           if max(max(actualImg)) == 1
              actualImg = uint8(actualImg*255);
           end
           oldOriginalImg = actualImg;
           actualImg = zeros(size(actualImg,1), size(actualImg,2), 3, 'uint8');
           for bandItr = 1:3
              actualImg(:,:,bandItr) = oldOriginalImg; 
           end
        end
        
        reconstructedMask = zeros(size(img,1), size(img,2), filtBandCount, 'uint8');
        labeledReconstructedMask = zeros(size(img,1), size(img,2));
        sizeOfImage = [size(img,1), size(img,2)];
        
        %% Go over each leaf in the graph and write its mask to the reconstructed image.
        nodes = imageNodeSets{fileItr};
        
        % If no nodes exist, do nothing.
        if isempty(nodes)
           continue; 
        end
        
        % Learn if this image is a background image.
        isBackground = ~nodes(1).sign;
        
        if isBackground
            outputDir = [backgroundDir '/' fileName];
        else 
            outputDir = [outputTempDir '/' fileName];
        end
        if ~exist(outputDir, 'dir')
           mkdir(outputDir); 
        end
        
        % Get the correct node set to reconstruct.
        if strcmp(reconstructionType, 'true')
            nodeReconInfo = [cat(1,nodes.labelId), cat(1,nodes.position)];
        else
            nodeReconInfo = leafNodes(:,1:3); %#ok<PFBNS>
        end
        
        %% Reconstruct each node.
        for nodeItr = 1:numel(nodes)
            % Fetch the node list to reconstruct.
            if strcmp(reconstructionType, 'true')
                reconstructedNodes = nodeItr;
            else
                reconstructedNodes = nodes(nodeItr).leafNodes;
            end
            
            %% Process each reconstructed node, and write them to a mask if necessary.
            reconstructedNodes = nodeReconInfo(reconstructedNodes,:);
            for reconNodeItr = 1:size(reconstructedNodes,1)
                % Read the mask here, and crop it if necessary.
                nodeMask = vocabMasks{reconstructedNodes(reconNodeItr,1)}; %#ok<PFBNS>
                position = reconstructedNodes(reconNodeItr,2:3);
                
                % Learn printed dimensions.
                halfSize = ceil((size(nodeMask)-1)/2);
                halfSize = halfSize(1:2);
                otherHalfSize = ([size(nodeMask,1), size(nodeMask,2)] - halfSize) -1; 
                imageSize = size(reconstructedMask);

                %% If patch is out of bounds, do nothing.
                if imageSize(1) < position(1)+halfSize(1) || ...
                        1 > (position(1)-halfSize(1)) || ...
                        imageSize(2) < (position(2)+halfSize(2)) ||...
                        1 > (position(2)-halfSize(2))
                    continue;
                end
                
                %% Write to reconstruction mask.
                 writtenMask = nodeMask;
                 reconstructedMask((position(1)-halfSize(1)):(position(1)+otherHalfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+otherHalfSize(2)),:) = ... 
                     max(reconstructedMask((position(1)-halfSize(1)):(position(1)+otherHalfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+otherHalfSize(2)),:), ...
                     writtenMask);
     
                 % First, write the node label to the labeled mask.
                 reconstructedPatch = labeledReconstructedMask((position(1)-halfSize(1)):(position(1)+otherHalfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+otherHalfSize(2)));
                 
                 avgNodeMask = mean(nodeMask,3);
                 reconstructedPatch(avgNodeMask > 10) = nodes(nodeItr).labelId;
                 labeledReconstructedMask((position(1)-halfSize(1)):(position(1)+otherHalfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+otherHalfSize(2))) = reconstructedPatch;
            end
            
            %% Mark the center of each realization with its label id.
            labeledReconstructedMask((nodes(nodeItr).position(1)-1):(nodes(nodeItr).position(1)+1), ...
                (nodes(nodeItr).position(2)-1):(nodes(nodeItr).position(2)+1)) = nodes(nodeItr).labelId;
        end
        
        %% Write the reconstructed mask to the output.
        % Add some random colors to make each composition look different, 
        % and overlay the gabors with the original image.
        rgbImg = label2rgb(labeledReconstructedMask, 'jet', 'k');
        %% Write the original image to a mask.
        if size(reconstructedMask,3)>1
            assignedBands = 1:size(reconstructedMask,3);
        else
            assignedBands = ones(size(rgbImg,3),1);
        end
        for bandItr = 1:size(rgbImg,3)
            rgbImg(:,:,bandItr) = uint8(double(rgbImg(:,:,bandItr)) .* (double(reconstructedMask(:,:,assignedBands(bandItr)))/255)) + ...
            actualImg(:,:,bandItr) .* uint8(reconstructedMask(:,:,assignedBands(bandItr))==0);
        end 

        %% Add edges to the image for visualization.
        edgeImg = zeros(sizeOfImage);
        edgeRgbImg = rgbImg;
        for nodeItr = 1:numel(nodes)
            edgeImg((nodes(nodeItr).position(1)-2):(nodes(nodeItr).position(1)+2), ...
                (nodes(nodeItr).position(2)-2):(nodes(nodeItr).position(2)+2)) = nodes(nodeItr).labelId;
            edges = nodes(nodeItr).adjInfo;
            if isempty(edges)
               continue;
            end
            edges = [edges(:,1:2) - nodeOffset, edges(:,3)];
            if ~isempty(edges)
                for edgeItr = 1:size(edges,1)
                   edgeIdx = drawline(double(nodes(edges(edgeItr,1)).position), double(nodes(edges(edgeItr,2)).position), sizeOfImage);
                   edgeImg(edgeIdx) = edges(edgeItr,3);
                end
            end
        end
        edgeImg = label2rgb(edgeImg, 'jet', 'k', 'shuffle');
        edgeRgbImg = max(edgeRgbImg, edgeImg);

        if levelItr>1
            imwrite(rgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) '_' reconstructionType '.png']);
            imwrite(edgeImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'onlyEdges.png']);
            imwrite(reconstructedMask, [outputDir, '/' fileName '_level' num2str(levelItr) 'clean.png']);
            imwrite(edgeRgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'edges_' reconstructionType '.png']);
        else
            imwrite(rgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) '_' reconstructionType '.png']);
            imwrite(edgeImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'onlyEdges.png']);
            imwrite(edgeRgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'edges_' reconstructionType  '.png']);
            imwrite(reconstructedMask, [outputDir, '/' fileName '_level' num2str(levelItr) 'clean.png']);
        end
    end
end