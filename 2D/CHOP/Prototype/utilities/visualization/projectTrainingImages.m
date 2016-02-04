%> Name: projectImages
%>
%> Description: Given the inference of parts in a dataset, this function
%> backprojects detected object models in real images. The backprojection
%> starts from the top level detections in every image, and continues all
%> the way to the bottom level. The same operation is performed both
%> training and test images.
%>
%> @param datasetName Name of the dataset to work on. 
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.11.2015
function [ ] = projectTrainingImages( fileList, vocabulary, mainGraph, levelItr, options)
    graphLevel = mainGraph{levelItr};
    imageIds = [graphLevel.imageId];
    realLabelIds = double(cat(1, graphLevel.realLabelId));
    precisePositions = double(cat(1, graphLevel.precisePosition));
    fastFlag = true;
    minIndividualPrint = 4;
    
    parfor imgItr = 1:max(imageIds)
       % Obtain image-specific realizations.
       idx = imageIds == imgItr;
       
       % If no nodes are present, move on.
       if nnz(idx) == 0
            continue;
       end
       imagePrecisePositions = precisePositions(idx,:);

       % Create a folder if needed.
       [~, fileName, ~] = fileparts(fileList{imgItr});
       orgImg = imread(fileList{imgItr});
       imgFolder = [pwd '/output/' options.datasetName '/reconstruction/train/' fileName];

       % Backproject from all possible levels.
       curExportArr = [realLabelIds(idx,:), imagePrecisePositions, repmat(levelItr, nnz(idx),1)];

       % Now, we get the top realizations and backproject to the original
       % image.
       
       level1Nodes = projectNode(curExportArr, vocabulary, 0, 'modal');

       % For visualization, overlay the original image with reconstructed nodes.
       [muImg, ~] = obtainPoE(level1Nodes, options.imageSize, options, fastFlag);
       muImg = muImg / max(max(max(muImg)));
       imwrite(muImg, [imgFolder '/' fileName '_level' num2str(levelItr) 'imagination.png']);
       imwrite(orgImg, [imgFolder '/' fileName '_original.png']);
       
       if levelItr >= minIndividualPrint
            imaginationFolder = [imgFolder '/' fileName '_level' num2str(levelItr) 'imaginations'];
            if ~exist(imaginationFolder, 'dir')
                 mkdir(imaginationFolder);
            end
            muImg = zeros(options.imageSize);
            varImg = zeros(options.imageSize);
            for nodeItr = 1:size(curExportArr,1)
                 partImg = imread([options.debugFolder '/level' num2str(levelItr) '/modalProjection/' num2str(curExportArr(nodeItr,1)) '.png']);
                 level1Nodes = projectNode(curExportArr(nodeItr,:), vocabulary, 0, 'modal');
                 [muImg, varImg] = obtainPoE(level1Nodes, options.imageSize, options, fastFlag, muImg, varImg);
                 partImg = double(partImg);
                 partImg = partImg / max(max(partImg));
                 
                 % If part image is too big, downsample it.
                 firstRatio = (size(muImg,1) / 5) / size(partImg,1);
                 secRatio = (size(muImg,2) / 5) / size(partImg,2);
                 if firstRatio < 1 || secRatio < 1
                      partImg = imresize(partImg, min(firstRatio, secRatio));
                 end
                 
                 % Print part image to the matrix.
                 visImg = muImg / max(max(max(muImg)));
                 visImg(1:size(partImg,1)+2, 1:size(partImg,2)+2, :) = 1;
                 visImg(2:size(partImg,1)+1, 2:size(partImg,2)+1, :) = partImg;
                 imwrite(visImg, [imgFolder '/' fileName '_level' num2str(levelItr) 'imaginations/realization_' num2str(nodeItr) '.png']);
            end
       end
    end