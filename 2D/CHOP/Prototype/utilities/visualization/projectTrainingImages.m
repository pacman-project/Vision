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
%    minIndividualPrint = 4;
    batchFlag = true;
    
    if strcmp(options.filterType, 'gabor')
        inhibitionHalfSize = options.gabor.inhibitionRadius;
        stride = options.gabor.stride;
    else
        inhibitionHalfSize = options.auto.inhibitionRadius;
        stride = options.auto.stride;
    end
    
    
    minPixelValue = 1/255;
    % Create filters.
    filters = options.filters;
    filters = cellfun(@(x) (x - min(min(x))) / (max(max(x)) - min(min(x))), filters, 'UniformOutput', false);
    for filterItr = 1:numel(filters)
         curFilter = filters{filterItr};
         curFilter(curFilter<minPixelValue) = minPixelValue;
         filters{filterItr} = curFilter;
    end
    visFilters = filters;
    
    % As a case study, we replace gabors with 1D gaussian filters
    % stretched.
    filterSize = size(filters{1},1) + 16;
    vals = normpdf(1:filterSize, (filterSize+1)/2, 5);
    vals = vals/max(vals);
    firstFilter = repmat(vals, filterSize, 1);
    angle = 180/numel(filters);
    for filterItr = 0:(numel(filters)-1)
         curAngle = -angle * filterItr;
         curFilter = imrotate(firstFilter, curAngle, 'bilinear', 'crop');
         curFilter(curFilter<minPixelValue) = minPixelValue;
         filters{filterItr+1} = curFilter;
    end
    
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
       prevRFSize = options.receptiveFieldSize * stride * (options.poolDim ^ (levelItr-3)) * (inhibitionHalfSize+1);
       [muImg, ~] = optimizeImagination(curExportArr, vocabulary, options.imageSize, prevRFSize, filters, visFilters, 1, batchFlag, options.datasetName, fileName);
%       level1Nodes = projectNode(curExportArr, vocabulary, 0, 'modal');

       % For visualization, overlay the original image with reconstructed nodes.
%       [muImg, ~] = obtainPoE(level1Nodes, [], [], [], options.imageSize, filters);
%       muImg = muImg / max(max(max(muImg)));
       imwrite(muImg, [imgFolder '/' fileName '_level' num2str(levelItr) 'optimizedImagination.png']);
       imwrite(orgImg, [imgFolder '/' fileName '_original.png']);
       
%        if levelItr >= minIndividualPrint
%             imaginationFolder = [imgFolder '/' fileName '_level' num2str(levelItr) 'imaginations'];
%             if ~exist(imaginationFolder, 'dir')
%                  mkdir(imaginationFolder);
%             end
%             muImg = zeros(options.imageSize);
%             varImg = zeros(options.imageSize);
%             likelihoodImg = zeros(options.imageSize);
%             for nodeItr = 1:size(curExportArr,1)
%                  partImg = imread([options.debugFolder '/level' num2str(levelItr) '/modalProjection/' num2str(curExportArr(nodeItr,1)) '.png']);
%                  level1Nodes = projectNode(curExportArr(nodeItr,:), vocabulary, 0, 'modal');
%                  
%                  [muImg, varImg, likelihoodImg, ~] = obtainPoE(level1Nodes, muImg, varImg, likelihoodImg, options.imageSize, filters);
%                  partImg = double(partImg);
%                  partImg = partImg / max(max(partImg));
%                  
%                  % If part image is too big, downsample it.
%                  firstRatio = (size(muImg,1) / 5) / size(partImg,1);
%                  secRatio = (size(muImg,2) / 5) / size(partImg,2);
%                  if firstRatio < 1 || secRatio < 1
%                       partImg = imresize(partImg, min(firstRatio, secRatio));
%                  end
%                  
%                  % Print part image to the matrix.
%                  visImg = muImg / max(max(max(muImg)));
%                  visImg(1:size(partImg,1)+2, 1:size(partImg,2)+2, :) = 1;
%                  visImg(2:size(partImg,1)+1, 2:size(partImg,2)+1, :) = partImg;
%                  imwrite(visImg, [imgFolder '/' fileName '_level' num2str(levelItr) 'imaginations/realization_' num2str(nodeItr) '.png']);
%             end
%        end
    end
    
    % Finally, create output images.
    for imgItr = 1:max(imageIds)
         [~, fileName, ~] = fileparts(fileList{imgItr});
         folderName = [pwd '/optimization/' options.datasetName '/level' num2str(levelItr) '/' fileName '_sample' num2str(1)];
         load([folderName '.mat'], 'imageArr', 'posLikelihoodArr', 'poeLikelihoodArr', 'diffImageArr');
         posPadding = (max(posLikelihoodArr) - min(posLikelihoodArr))/4;
         posLimits = [min(posLikelihoodArr) - posPadding, max(posLikelihoodArr) + posPadding];
         if numel(unique(posLimits)) == 1
              posLimits = [posLimits(1)-1, posLimits(1)+1];
         end
         poePadding = (max(poeLikelihoodArr) - min(poeLikelihoodArr))/4;
         poeLimits = [min(poeLikelihoodArr) - poePadding, max(poeLikelihoodArr) + poePadding];
         if numel(unique(poeLimits)) == 1
              poeLimits = [(poeLimits(1)-1), (poeLimits(1) +1)];
         end
         fileName = [folderName '.gif'];
         for stepItr = 1:numel(imageArr)
              figure('Visible', 'off'), hold on;
              axis square
              subplot(2,2,1), imshow(imageArr{stepItr});
              title('Imagined Data')

              subplot(2,2,2)
              title('Likelihood differences')
              if ~isempty(diffImageArr{stepItr})
                   imshow(diffImageArr{stepItr});
              end

              subplot(2,2,3), plot(1:stepItr, posLikelihoodArr(1:stepItr));
              ylim(posLimits);
              if stepItr>1
                    xlim([1, numel(imageArr)]);
              end
              title('Change in position likelihood')

              subplot(2,2,4), plot(1:stepItr, poeLikelihoodArr(1:stepItr));
              ylim(poeLimits);
              if stepItr>1
                    xlim([1, numel(imageArr)]);
              end
              title('Product of experts likelihood')

              hold off;
              saveas(gcf,  [folderName '_temp.png']);
              im=imread([folderName '_temp.png']);
              [imind,cm] = rgb2ind(im, 256);

              if stepItr == 1
                   imwrite(imind, cm, fileName, 'LoopCount', inf, 'DelayTime',2);
              else
                   imwrite(imind, cm, fileName, 'WriteMode', 'append');
              end
              close(gcf);
         end
    end
    
    