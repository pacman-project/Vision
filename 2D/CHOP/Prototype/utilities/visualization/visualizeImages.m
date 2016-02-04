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
function [ ] = visualizeImages( fileList, vocabLevel, graphLevel, representativeNodes, allNodeInstances, leafNodes, leafNodeCoords, levelItr, options, type )
    outputTempDir = [options.outputFolder '/reconstruction/' type];
    backgroundDir = [options.outputFolder '/reconstruction/' type '/' options.backgroundClass];
    processedFolder = options.processedFolder;   
    reconstructionType = options.reconstructionType;
    instancePerNode = options.vis.instancePerNode-1;
    printTrainRealizations = options.vis.printTrainRealizations;
    numberOfVocabLevelNodes = double(max([vocabLevel.label]));
    activations = cat(1, graphLevel.activation);
    
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
            vocabMasks{fileItr} = imread([options.currentFolder '/debug/' options.datasetName ...
                '/level1/reconstruction/' num2str(fileItr) '.png']);
        end
    else
        vocabMaskList = fuf([options.currentFolder '/debug/' options.datasetName '/level' num2str(usedLevel) '/reconstruction/*.png'], 1, 'detail');
        minusList = fuf([options.currentFolder '/debug/' options.datasetName '/level' num2str(usedLevel) '/reconstruction/*_comp.png'], 1, 'detail');
        minusList = [minusList; fuf([options.currentFolder '/debug/' options.datasetName '/level' num2str(usedLevel) '/reconstruction/*_uni.png'], 1, 'detail')];
        numberOfMasks = numel(vocabMaskList) - numel(minusList);
        for fileItr = 1:numberOfMasks
            vocabMasks{fileItr} = imread([options.currentFolder '/debug/' options.datasetName ...
                '/level' num2str(usedLevel) '/reconstruction/' num2str(fileItr) '.png']);
        end
    end
    
    %% Create folders to put the cropped original images for each realization in.
    usedChildren = ones(numel(vocabLevel), 1) * instancePerNode;
    croppedOrgFolder = [options.currentFolder '/debug/' options.datasetName '/level' num2str(levelItr) '/cropped'];
    if ~exist(croppedOrgFolder, 'dir')
       mkdir(croppedOrgFolder); 
    end
    
    %% Put realizations into distinct sets so that each image has its own nodes in a set.
    imageIds = [graphLevel.imageId]';
    imageNodeSets = cell(numel(fileList),1);
    imageNodeIds = cell(numel(fileList), 1);
    for fileItr = 1:numel(fileList)
        imageNodeSets(fileItr) = {graphLevel(imageIds == fileItr)};
        imageNodeIds(fileItr) = {find(imageIds == fileItr)};
    end
    
    %% Go over the list of images and run reconstruction.
    for fileItr = 1:numel(fileList)
        nodeOffset = numel(find(imageIds<fileItr));
        relevantNodeIds = imageNodeIds{fileItr};
        
        %% Learn the size of the original image, and allocate space for new mask.
        [~, fileName, ~] = fileparts(fileList{fileItr});
        img = imread([processedFolder '/' fileName '.png']);
        actualImg = img;
        
        % If this is a depth image, we need to reduce the data
        % to 8 bit.
        if isa(actualImg, 'uint16')
             actualImg = uint8(double(actualImg)/256);
        end
        
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
        
        reconstructedMask = zeros(size(img,1), size(img,2), filtBandCount, 'double');
        reconstructedMaskCounts = zeros(size(img,1), size(img,2), 'double');
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
            nodeReconInfo = [cat(1,nodes.labelId), cat(1,nodes.precisePosition)];
        else
            nodeReconInfo = [leafNodes(:,1), leafNodeCoords]; %#ok<PFBNS>
        end
        
        %% Reconstruct each node.
        for nodeItr = 1:numel(nodes)
            activation = nodes(nodeItr).activation;
            
            % Fetch the node list to reconstruct.
            if strcmp(reconstructionType, 'true')
                reconstructedNodes = nodeItr;
            else
                reconstructedNodes = nodes(nodeItr).leafNodes;
            end
            
            %% Process each reconstructed node, and write them to a mask if necessary.
            reconstructedNodes = nodeReconInfo(reconstructedNodes,:);
            minX = Inf;
            maxX = -1;
            minY = Inf;
            maxY = -1;
            for reconNodeItr = 1:size(reconstructedNodes,1)
                % Read the mask here, and crop it if necessary.
                nodeMask = vocabMasks{reconstructedNodes(reconNodeItr,1)}; %#ok<PFBNS>
                position = reconstructedNodes(reconNodeItr,2:3);
                
                % Learn printed dimensions.
                halfSize = floor((size(nodeMask)-1)/2);
                halfSize = halfSize(1:2);
                otherHalfSize = ([size(nodeMask,1), size(nodeMask,2)] - halfSize) - 1; 
                imageSize = size(reconstructedMask);

                %% If patch is out of bounds, do nothing.
                if imageSize(1) < position(1)+halfSize(1) || ...
                        1 > (position(1)-halfSize(1)) || ...
                        imageSize(2) < (position(2)+halfSize(2)) ||...
                        1 > (position(2)-halfSize(2))
                    continue;
                end
                
                %% Write to reconstruction mask.
                 if minX > (position(1)-halfSize(1))
                     minX = (position(1)-halfSize(1));
                 end
                 if maxX < (position(1)+otherHalfSize(1))
                     maxX = (position(1)+otherHalfSize(1));
                 end
                 if minY > (position(2)-halfSize(2))
                     minY = (position(2)-halfSize(2));
                 end
                 if maxY < (position(2)+otherHalfSize(2))
                     maxY = (position(2)+otherHalfSize(2));
                 end
                 writtenMask = double(nodeMask);
                 writtenMaskMap = double(mean(writtenMask,3) > 0);
                 reconstructedMask((position(1)-halfSize(1)):(position(1)+otherHalfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+otherHalfSize(2)),:) = ... 
                     reconstructedMask((position(1)-halfSize(1)):(position(1)+otherHalfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+otherHalfSize(2)),:) + ...
                     writtenMask;
                 reconstructedMaskCounts((position(1)-halfSize(1)):(position(1)+otherHalfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+otherHalfSize(2)),:) = ...
                     reconstructedMaskCounts((position(1)-halfSize(1)):(position(1)+otherHalfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+otherHalfSize(2)),:) + ...
                     writtenMaskMap;
                     
                 % First, write the node label to the labeled mask.
                 reconstructedPatch = labeledReconstructedMask((position(1)-halfSize(1)):(position(1)+otherHalfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+otherHalfSize(2)));
                 
                 avgNodeMask = mean(nodeMask,3);
                 reconstructedPatch(avgNodeMask > 10) = nodes(nodeItr).labelId;
                 labeledReconstructedMask((position(1)-halfSize(1)):(position(1)+otherHalfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+otherHalfSize(2))) = reconstructedPatch;
            end
            
            %% Write original image's cropped area to a file.
            if ~isempty(allNodeInstances)
                 representativeNode = representativeNodes(nodes(nodeItr).labelId);
                 nodeInstances = allNodeInstances{representativeNode};
                 if ismember(relevantNodeIds(nodeItr), nodeInstances(:,1))
                     imgId = find(nodeInstances(:,1) == relevantNodeIds(nodeItr)) - 1;
                     imwrite(actualImg(nodeInstances(imgId+1,2):nodeInstances(imgId+1,4), ...
                         nodeInstances(imgId+1,3):nodeInstances(imgId+1,5), :), [croppedOrgFolder '/' num2str(representativeNode) '_' num2str(imgId) '.png']);
                    usedChildren(representativeNode) = usedChildren(representativeNode) -1;
                 end
                 
                 nodeInstances = allNodeInstances{nodes(nodeItr).realLabelId};
                 if ismember(relevantNodeIds(nodeItr), nodeInstances(:,1)) && nodes(nodeItr).realLabelId ~= representativeNode
                     imgId = find(nodeInstances(:,1) == relevantNodeIds(nodeItr)) - 1;
                     imwrite(actualImg(nodeInstances(imgId+1,2):nodeInstances(imgId+1,4), ...
                         nodeInstances(imgId+1,3):nodeInstances(imgId+1,5), :), [croppedOrgFolder '/' num2str(nodes(nodeItr).realLabelId) '_' num2str(imgId) '.png']);
                    usedChildren(nodes(nodeItr).realLabelId) = usedChildren(nodes(nodeItr).realLabelId) -1;
                 end
            end
            
            %% Mark the center of each realization with its label id.
            labeledReconstructedMask(round(nodes(nodeItr).precisePosition(1)-1):round(nodes(nodeItr).precisePosition(1)+1), ...
                round(nodes(nodeItr).precisePosition(2)-1):round(nodes(nodeItr).precisePosition(2)+1)) = nodes(nodeItr).labelId;
        end
       
        % Normalize and obtain final reconstructed image.
        for bandItr = 1:size(reconstructedMask,3)
             relevantImg = reconstructedMask(:,:,bandItr);
             relevantImg(reconstructedMaskCounts > 0) = relevantImg(reconstructedMaskCounts > 0) ./ ...
                 reconstructedMaskCounts(reconstructedMaskCounts > 0);
             reconstructedMask(:,:,bandItr) = relevantImg;
        end
        reconstructedMask = uint8(round(reconstructedMask));
        
        if strcmp(type, 'test') || (strcmp(type, 'train') && printTrainRealizations)
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
            
            nodeImg = zeros(sizeOfImage);
            centerSize = 1;
            for nodeItr = 1:numel(nodes)
                edgeImg(round(nodes(nodeItr).precisePosition(1)-centerSize):round(nodes(nodeItr).precisePosition(1)+centerSize), ...
                    round(nodes(nodeItr).precisePosition(2)-centerSize):round(nodes(nodeItr).precisePosition(2)+centerSize)) = nodes(nodeItr).labelId;
                nodeImg(round(nodes(nodeItr).precisePosition(1)-centerSize):round(nodes(nodeItr).precisePosition(1)+centerSize), ...
                    round(nodes(nodeItr).precisePosition(2)-centerSize):round(nodes(nodeItr).precisePosition(2)+centerSize)) = nodes(nodeItr).labelId;
            end
            
            for nodeItr = 1:numel(nodes)
                edges = nodes(nodeItr).adjInfo;
                if ~isempty(edges)
                    edges = [edges(:,1:2) - nodeOffset, edges(:,3)];
                    if ~isempty(edges)
                        for edgeItr = 1:size(edges,1)
                           edgeIdx = drawline(double(nodes(edges(edgeItr,1)).precisePosition), double(nodes(edges(edgeItr,2)).precisePosition), sizeOfImage);
                           edgeImg(edgeIdx) = edges(edgeItr,3);
                        end
                    end
                end
            end
            
            % 
            edgeImg = label2rgb(edgeImg, 'jet', 'k', 'shuffle');
            nodeImg = label2rgb(nodeImg, 'jet', 'k', 'shuffle');
            edgeRgbImg = max(edgeRgbImg, edgeImg);

            % 
            if levelItr>1
                imwrite(rgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) '_' reconstructionType '.png']);
                imwrite(edgeImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'onlyEdges.png']);
                imwrite(reconstructedMask, [outputDir, '/' fileName '_level' num2str(levelItr) 'projection.png']);
                imwrite(nodeImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'realizations.png']);
                imwrite(edgeRgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'edges_' reconstructionType '.png']);
            else
               imwrite(rgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) '_' reconstructionType '.png']);
                imwrite(edgeImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'onlyEdges.png']);
               imwrite(edgeRgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'edges_' reconstructionType  '.png']);
                imwrite(reconstructedMask, [outputDir, '/' fileName '_level' num2str(levelItr) 'projection.png']);
                imwrite(nodeImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'realizations.png']);
            end
        end
    end
end
