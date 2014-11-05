%> Name: visualizeLevel
%>
%> Description: Given the current and previous vocabulary level, this
%> function visualizes the current vocabulary level as a separate patche for
%> each word in the vocabulary level.
%>
%> @param currentLevel Current vocabulary level.
%> @param graphLevel Current graph level.
%> @param levelId Identifier of the current level.
%> @param numberOfPrevNodes Number of nodes in previous vocabulary level.
%> @param options Program options.
%> @param isRedundant If currentLevel consists of redundant compositions,
%> set 1. Otherwise set 0.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 10.02.2014
%> Redundant vocabulary output option added. 10.05.2014
function [] = visualizeCroppedImgs( currentLevel, levelId, options)
    % Read options to use in this file.
    currentFolder = options.currentFolder;
    datasetName = options.datasetName;
    instancePerNode = options.vis.instancePerNode;
    instanceImgDim = round(sqrt(instancePerNode));
    
    if isempty(currentLevel)
       return; 
    end
    
    %% Create output folder structure.
    croppedDir = [currentFolder '/debug/' datasetName '/level' num2str(levelId) '/cropped/'];
    
    %% Combine all compositions and show them within a single image.
    % Learn number of rows/columns.
    numberOfNodes = numel(currentLevel);
    colImgCount = ceil(sqrt(numberOfNodes));
    rowImgCount = ceil(numberOfNodes/colImgCount);

    % Read cropped images.
    nodeImgs = cell(numberOfNodes,1);
    for nodeItr = 1:numel(currentLevel)
        instanceImgs = cell(instancePerNode-1, 1);
        for instItr = 1:(instancePerNode-1)
            filePath = [croppedDir '/' num2str(nodeItr) '_' num2str(instItr) '.png'];
            if exist([croppedDir '/' num2str(nodeItr) '_' num2str(instItr) '.png'], 'file')
                img = imread(filePath);
                instanceImgs(instItr) = {img};
            else
                break;
            end
        end
        instanceImgs = instanceImgs(cellfun(@(x) ~isempty(x), instanceImgs));
        nodeImgs(nodeItr) = {instanceImgs};
    end
    
    % Learn the right mask sizes and band counts.
    compMaskSize = [1, 1];
    dim3 = 1;
    for nodeItr = 1:numberOfNodes
        instanceImgs = nodeImgs{nodeItr};
        if ~isempty(instanceImgs)
            dim3 = size(instanceImgs{1},3);
            instanceImgSizes = cellfun(@(x) [size(x,1), size(x,2)], instanceImgs, 'UniformOutput', false);
            compMaskSize = max(compMaskSize, max(cat(1, instanceImgSizes{:})));
        end
    end
    
    % Put the part to the instance list, too.
    for nodeItr = 1:numberOfNodes
        instanceImgs = nodeImgs{nodeItr};
        filterImg = imread([currentFolder '/debug/' datasetName '/level' num2str(levelId) '/reconstruction/' num2str(nodeItr) '.png']);
        newSize=min([size(filterImg,1), size(filterImg,2)], compMaskSize);
        filterImg = filterImg(1:newSize(1), 1:newSize(2), :);
        if size(filterImg,3) < dim3
            filterMultipleImg = zeros(size(filterImg,1), size(filterImg,2), dim3, 'uint8');
            for dimItr = 1:dim3
                filterMultipleImg(:,:,dimItr) = filterImg;
            end
        end
        instanceImgs = [{filterMultipleImg}; instanceImgs];
        nodeImgs(nodeItr) = {instanceImgs};
    end
    
    % Make mask sizes uniform and write them all back.
    for nodeItr = 1:numberOfNodes
        instanceImgs = nodeImgs{nodeItr};
        for instItr = 1:numel(instanceImgs)
            tempMask2 = instanceImgs{instItr};
            finalTempMask = zeros([compMaskSize, size(tempMask2,3)], 'uint8');
            margins = (compMaskSize - [size(tempMask2,1), size(tempMask2,2)])/2;
            finalTempMask((floor(margins(1))+1):(end-ceil(margins(1))), ...
                (floor(margins(2))+1):(end-ceil(margins(2))), :) = tempMask2;
            % A make-up to fill in NaNs (empty points).
            fillInValue = median(double(finalTempMask(finalTempMask>0 & finalTempMask<255)));
            finalTempMask(finalTempMask == 0) = fillInValue;
            instanceImgs{instItr} = finalTempMask;
            imwrite(finalTempMask, [croppedDir '/' num2str(nodeItr) '_' num2str(instItr) '.png']);
        end
        nodeImgs{nodeItr} = instanceImgs;
    end

    % Using the maximum dimensions, transform each composition image to the
    % same size. 
    overallInstanceImage = NaN((rowImgCount * instanceImgDim)*(compMaskSize(1)+1)+1, colImgCount * instanceImgDim * (compMaskSize(2)+1)+1, dim3);
 %   overallInstanceRealImage = NaN((rowImgCount * instanceImgDim)*(compMaskSize(1)+1)+1, colImgCount * instanceImgDim * (compMaskSize(2)+1)+1, dim3);
    for nodeItr = 1:numberOfNodes
        instanceImgs = nodeImgs{nodeItr};
        for instItr = 1:numel(instanceImgs)
            compFinalMask = instanceImgs{instItr};

            % Add the composition's mask to the overall mask image.
            rowStart2= 2 + floor((nodeItr-1)/colImgCount)*(compMaskSize(1)+1) *instanceImgDim ;
            colStart2 = 2 + rem(nodeItr-1, colImgCount) * (compMaskSize(2)+1) * instanceImgDim;
            %We're writing sample instances to other images. Find where to
            %write them and put them in their location.
            rowInstStart = floor((instItr - 1)/instanceImgDim)*(compMaskSize(1)+1);
            colInstStart = rem(instItr-1, instanceImgDim) * (compMaskSize(2)+1);
            overallInstanceImage((rowStart2 + rowInstStart):((rowStart2 + rowInstStart)+compMaskSize(1)-1), ...
                (colStart2+colInstStart):((colStart2+colInstStart)+compMaskSize(2)-1), :) = compFinalMask;
        end
  %      imwrite(instanceImgs{1}, [reconstructionDir num2str(nodeItr) '_uni.png']);
    end

    clear instanceImgs nodeImgs setImgs;

    % A final make up in order to separate masks from each other by 1s.
    whiteRowIdx = 1:(compMaskSize(1)+1) *instanceImgDim:size(overallInstanceImage,1);
    whiteColIdx = 1:(compMaskSize(2)+1) *instanceImgDim:size(overallInstanceImage,2);
    overallInstanceImage(isnan(overallInstanceImage)) = 0;
    overallInstanceImage(:, whiteColIdx, :) = 255;
    overallInstanceImage(whiteRowIdx, :, :) = 255;
    overallInstanceImage = uint8(overallInstanceImage);

    % Then, write the compositions the final image.
    imwrite(overallInstanceImage, [currentFolder '/debug/' datasetName '/level' num2str(levelId) '_vb_cropped.png']);
end



