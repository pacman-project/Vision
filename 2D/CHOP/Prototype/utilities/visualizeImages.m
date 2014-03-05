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
function [ ] = visualizeImages( fileList, vocabLevel, graphLevel, leafNodes, levelItr, options, type )
    outputDir = [options.outputFolder '/reconstruction/' type];
    originalDir = [options.outputFolder '/original/'];
    if ~exist(outputDir, 'dir')
       mkdir(outputDir); 
    end
    processedDir = options.processedFolder;
    reconstructionType = options.reconstructionType;
    imageReconstructionType = options.imageReconstructionType;
    
    if strcmp(imageReconstructionType, 'individual') && levelItr < options.minIndividualReconstructionLevel
        return;
    end
    
    
    %% Depending on the reconstruction type, we read masks to put on correct positions in images.
    if strcmp(options.reconstructionType, 'true')
        % Read masks of compositions in current level.
        usedLevel = levelItr;
    else
        % Read level 1 into separate masks.
        usedLevel = 1;
    end
    
    vocabMaskList = fuf([options.currentFolder '/debug/' options.datasetName '/level' num2str(usedLevel) '/reconstruction/*.png'], 1, 'detail');
    numberOfMasks = numel(vocabMaskList);
    vocabMasks = cell(numberOfMasks,1);

    for fileItr = 1:numberOfMasks
        vocabMasks{fileItr} = imread([options.currentFolder '/debug/' options.datasetName ...
            '/level' num2str(usedLevel) '/reconstruction/' num2str(fileItr) '.png']);
    end
    
    imageIds = [graphLevel.imageId]';
    imageNodeSets = cell(numel(fileList),1);
    for fileItr = 1:numel(fileList)
        imageNodeSets(fileItr) = {graphLevel(imageIds == fileItr)};
    end
    
    %% Go over the list of images and run reconstruction.
    for fileItr = 1:numel(fileList)
        nodeOffset = numel(find(imageIds<fileItr));
        %% Learn the size of the original image, and allocate space for new mask.
        [~, fileName, ~] = fileparts(fileList{fileItr});
        img = imread([processedDir '/' fileName '.png']);
        actualImg = imread([originalDir fileName '.png']);
        
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
            
        originalImg = actualImg;
        reconstructedMask = zeros(size(img,1), size(img,2), 'uint8');
        labeledReconstructedMask = zeros(size(img,1), size(img,2));
        sizeOfImage = [size(img,1), size(img,2)];
        
        %% Go over each leaf in the graph and write its mask to the reconstructed image.
        nodes = imageNodeSets{fileItr};
        
        % If no nodes exist, do nothing.
        if isempty(nodes)
           continue; 
        end
        
        % Get mdl scores 
        if levelItr == 1
            mdlScores = zeros(numel(vocabLevel),1);
            maxMdlScore = 1;
        else
            mdlScores = [vocabLevel.normMdlScore];
            maxMdlScore = max(mdlScores);
            if maxMdlScore == 0
                maxMdlScore = 1;
            end
        end
        
        
        % Get the correct node set to reconstruct.
        if strcmp(reconstructionType, 'true')
            nodeReconInfo = [{nodes.labelId}', {nodes.position}'];
        else
            nodeReconInfo = leafNodes(:,1:2);
        end
        
        %% Reconstruct each node.
        for nodeItr = 1:numel(nodes)
            % Fetch the node list to reconstruct.
            if strcmp(reconstructionType, 'true')
                reconstructedNodes = nodeItr;
            else
                reconstructedNodes = nodes(nodeItr).leafNodes;
            end
            
            % If each realization is to be written separately, raed a new
            % image each time.
            if strcmp(imageReconstructionType, 'individual')
                actualImg = originalImg;
                labeledReconstructedMask = zeros(size(img,1), size(img,2), 'uint8');
                reconstructedMask = labeledReconstructedMask;
                labeledReconstructedMask(end,end) = 255;
            end
            
            %% Process each reconstructed node, and write them to a mask if necessary.
            reconstructedNodes = nodeReconInfo(reconstructedNodes,:);
            for reconNodeItr = 1:size(reconstructedNodes,1)
                % Read the mask here, and crop it if necessary.
                if strcmp(imageReconstructionType, 'individual')
                    mdlScore = mdlScores(reconstructedNodes{reconNodeItr,1});
                end
                nodeMask = vocabMasks{reconstructedNodes{reconNodeItr,1}};
                position = reconstructedNodes{reconNodeItr,2};
                
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
                
                %% Write to reconstructin mask.
 %               if strcmp(imageReconstructionType, 'all')
                    reconstructedMask((position(1)-halfSize(1)):(position(1)+halfSize(1)), ...
                         (position(2)-halfSize(2)):(position(2)+halfSize(2)),:) = ... 
                         max(reconstructedMask((position(1)-halfSize(1)):(position(1)+halfSize(1)), ...
                         (position(2)-halfSize(2)):(position(2)+halfSize(2)),:), ...
                         nodeMask);
 %               end
                 % First, write the node label to the labeled mask.
                 reconstructedPatch = labeledReconstructedMask((position(1)-halfSize(1)):(position(1)+halfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+halfSize(2)));
                 
                 if strcmp(imageReconstructionType, 'all')
                     reconstructedPatch(nodeMask > 10) = nodeItr;
                 elseif levelItr>1
                     reconstructedPatch(nodeMask > 10) = round(255*(mdlScore/maxMdlScore));
                 end
                    
                 labeledReconstructedMask((position(1)-halfSize(1)):(position(1)+halfSize(1)), ...
                     (position(2)-halfSize(2)):(position(2)+halfSize(2))) = reconstructedPatch;
                 
            end
            %% Print this sub to a separate mask, if needed.
            if strcmp(imageReconstructionType, 'individual')
                rgbImg = label2rgb(labeledReconstructedMask, 'autumn', 'k', 'shuffle');

                for bandItr = 1:size(rgbImg,3)
                    rgbImg(:,:,bandItr) = uint8(double(rgbImg(:,:,bandItr)) .* (double(reconstructedMask))) + ...
                    actualImg(:,:,bandItr) .* uint8(reconstructedMask==0);
                end 

                % Write to output file.
                imwrite(rgbImg, ...
                    [outputDir, '/' fileName '_level' num2str(levelItr) '_realization_' num2str(nodeItr) ...
                    '_mdl_' num2str(mdlScore) '_' reconstructionType '.png']);
            end
        end
        
        if strcmp(imageReconstructionType, 'all')
            %% Write the reconstructed mask to the output.
            % Add some random colors to make each composition look different, 
            % and overlay the gabors with the original image.
            rgbImg = label2rgb(labeledReconstructedMask, 'jet', 'k', 'shuffle');
            %% Write the original image to a mask.
            for bandItr = 1:size(rgbImg,3)
                rgbImg(:,:,bandItr) = uint8(double(rgbImg(:,:,bandItr)) .* (double(reconstructedMask)/255)) + ...
                actualImg(:,:,bandItr) .* uint8(reconstructedMask==0);
            end 

            %% Add edges to the image for visualization.
            edgeImg = zeros(sizeOfImage);
            edgeRgbImg = rgbImg;
            for nodeItr = 1:numel(nodes)
                edges = nodes(nodeItr).adjInfo;
                if isempty(edges)
                   continue;
                end
                edges = edges(:,1:2) - nodeOffset;
                if ~isempty(edges)
                    for edgeItr = 1:size(edges,1)
                       edgeIdx = drawline(nodes(edges(edgeItr,1)).position, nodes(edges(edgeItr,2)).position, sizeOfImage);
                       edgeImg(edgeIdx) = 1;
                    end
                end
            end
            edgeImg = edgeImg > 0;

            for bandItr = 1:3
                if bandItr == 1
                    edgeRgbImg(:,:,bandItr) = max(rgbImg(:,:,bandItr), uint8(edgeImg * 255));
                else
                    edgeRgbImg(:,:,bandItr) = min(rgbImg(:,:,bandItr), uint8((~edgeImg) * 255));
                end
            end
        
            %% Print the images.
            if levelItr>1
                % Read first level to fill the voids left by missing filters.
                firstLevelMask = imread([outputDir, '/' fileName '_level1clean.png']);
                firstLevelImg = zeros(size(firstLevelMask,1), size(firstLevelMask,2), size(rgbImg,3), 'uint8');
                for bandItr = 1:size(rgbImg,3)
                    firstLevelImg(:,:,bandItr) = firstLevelMask;
                end

                % Combine both and write to output.
                imwrite(rgbImg + firstLevelImg, [outputDir, '/' fileName '_level' num2str(levelItr) '_' reconstructionType '.png']);
                imwrite(edgeImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'onlyEdges.png']);
                imwrite(edgeRgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'edges_' reconstructionType '.png']);
            else
                imwrite(rgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) '_' reconstructionType '.png']);
                imwrite(edgeImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'onlyEdges.png']);
                imwrite(edgeRgbImg, [outputDir, '/' fileName '_level' num2str(levelItr) 'edges_' reconstructionType  '.png']);
                imwrite(reconstructedMask, [outputDir, '/' fileName '_level' num2str(levelItr) 'clean.png']);
            end
        end
    end
end