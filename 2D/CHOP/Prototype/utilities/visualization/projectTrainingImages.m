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
    batchFlag = false;
    dummySigma = 0.15;
    angleStep = 3; % 22.5 for no rotations!
    
    if strcmp(options.filterType, 'gabor')
        inhibitionHalfSize = options.gabor.inhibitionRadius;
        stride = options.gabor.stride;
    else
        inhibitionHalfSize = options.auto.inhibitionRadius;
        stride = options.auto.stride;
    end
    
    minPixelValue = 1/255;
    
    % As a case study, we replace gabors with 1D gaussian filters
    % stretched.
    visFilterSize = size(options.filters{1},1);
    
    % Generate barrel filter (initial)
    filterSize = visFilterSize + 8;
    vals = normpdf(1:filterSize, (filterSize+1)/2, 3);
    vals = vals/max(vals);
    firstFilter = repmat(vals, filterSize, 1);
    
    % Obtain visual filter.
    firstVisFilter = options.filters{1};
    firstVisFilter = (firstVisFilter - min(min(firstVisFilter))) / (max(max(firstVisFilter)) - min(min(firstVisFilter)));
    
    % Generate rotated filters.
    numberOfOptFilters = round(180/angleStep);
    optFilters = cell(numberOfOptFilters,1);
    visFilters = cell(numberOfOptFilters,1);
    for filterItr = 0:((180-angleStep)/angleStep)
         curAngle = -angleStep * filterItr;
         curFilter = imrotate(firstFilter, curAngle, 'bilinear', 'crop');
         curVisFilter = imrotate(firstVisFilter, curAngle, 'bilinear', 'crop');
         curFilter(curFilter<minPixelValue) = minPixelValue;
         curVisFilter(curVisFilter<minPixelValue) = minPixelValue;
         curFilter = round(curFilter * 255)/255;
         curVisFilter = round(curVisFilter * 255)/255;
         optFilters{filterItr+1} = curFilter;
         visFilters{filterItr+1} = curVisFilter;
    end
        
    %% Now, we convert the filters to uint8.
    for filterItr = 1:numel(optFilters)
         visFilters{filterItr} = uint8(round(visFilters{filterItr} * 255));
         optFilters{filterItr} = uint8(round(optFilters{filterItr} * 255));
    end
    optFilters = cat(3, optFilters{:});
    visFilters = cat(3, visFilters{:});
    
    %% Here, we create a likelihood lookup table for predictions.
    likelihoodLookupTable = zeros(256);
    startVals = (-1/510:1/255:(1+1/510))';
    for itr = 0:255
         probs = normcdf(startVals, ones(257, 1) * (itr/255) , repmat(dummySigma,257,1));
 %        normProb = probs(end) - probs(1);
         normProb = 1;
         likelihoodLookupTable(itr+1,:) = (probs(2:end) - probs(1:(end-1))) / normProb;
    end
    likelihoodLookupTable = log(likelihoodLookupTable);

    %% Go through every image and find optimal version.
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
       [muImg, ~] = optimizeImagination(curExportArr, vocabulary, options.imageSize, prevRFSize, optFilters, visFilters, 1, batchFlag, options.datasetName, likelihoodLookupTable, fileName);

       % For visualization, overlay the original image with reconstructed nodes.
       imwrite(muImg, [imgFolder '/' fileName '_level' num2str(levelItr) 'optimizedImagination.png']);
       imwrite(orgImg, [imgFolder '/' fileName '_original.png']);
       
    end
end
    