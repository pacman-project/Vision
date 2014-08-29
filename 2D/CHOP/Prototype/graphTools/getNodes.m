%> Name: getInitialNodes
%>
%> Description: When given an image, this function is used to get initial
%> simple nodes by applying Gabor filters over edge mask of the image. Each
%> response peak is considered as a simple feature instance.
%>
%> @param img Input image
%> @param gtFileName An empty param either means the gt for that image is not
%> given, or it is of invalid size and cannot be used. If not empty, feel
%> free to use it in elimination of nodes (given the gt use is enabled in 
%> options).
%> @param options Program options.
%> 
%> @retval nodes The nodes to form further graphs.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.11.2013
%> Ver 1.1 on 03.12.2013 Response inhibition added.
%> Ver 1.2 on 12.01.2014 Comment changes to create unified code look.
%> Ver 1.3 on 17.01.2014 GT use implemented.
function [ nodes, smoothedImg ] = getNodes( img, gtFileName, options )
    %% Step 1: Get grayscaled image.
    if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
        if size(img,3)>1
            img = rgb2gray(img(:,:,1:3));
        end
    else
        filterMatrix = options.filterMatrix;
        deadFeatureStd = options.auto.deadFeatureStd;
        stride = options.auto.stride;
        whMat = options.auto.whMat;
        mu = options.auto.mu;
    end
    filterCount = numel(options.filters);
    filterSize = size(options.filters{1});
    filterBandSize = filterSize(1:2);
    
    %% Get gt info in the form of a mask.
    if options.useGT && ~isempty(gtFileName)
        gtMask = imread(gtFileName);
        if strcmp(options.gtType, 'contour')
            gtMask = imdilate(gtMask, strel('disk', options.contourGTNeighborhood, 8));
        else
            gtMask = imfill(gtMask, 'holes');
        end
    else
        gtMask = ones(size(img(:,:,1))) > 0;
    end
    
    % Apply denoising to get better responses.
    for bandItr = 1:size(img,3);
        myfilter = fspecial('gaussian',[3 3], 2);
        img(:,:,bandItr) = imfilter(img(:,:,bandItr), myfilter, 'replicate', 'same', 'conv');
        img(:,:,bandItr)=medfilt2(img(:,:,bandItr), [3,3]);
        img = double(img);
    end
    
    %% Get response by applying each filter to the image.
    responseImgs = zeros(size(img,1), size(img,2), filterCount);
    if strcmp(options.filterType, 'auto')
        dim1 = (size(img,1)-filterSize(1)+1);
        dim2 = (size(img,2)-filterSize(2)+1);

        % Implementing stride here. Some blocks will be skipped. 
        imgCols = zeros( ceil(dim1/stride) * ceil(dim2/stride), prod(filterSize));
        startIdx = 1;
        iterator = prod(filterBandSize)-1;
        
        % Get linear indices of correct columns, since not all blocks (each
        % corresponds a column) may be used, depending on the stride parameter.
        idx1 = 1:stride:dim1;
        idx2 = 1:stride:dim2;
        [p,q] = meshgrid(idx1, idx2);
        pairs = [p(:) q(:)];
        pairs = sortrows(pairs,2);
        validCols = sub2ind([dim1, dim2], pairs(:,1), pairs(:,2));
        
        for bandItr = 1:size(img,3)
            tempCols = im2col(img(:,:,bandItr), filterBandSize)';
            imgCols(:,startIdx:(startIdx+iterator)) = tempCols(validCols,:);
            startIdx = startIdx + iterator + 1;
        end
        muArr = repmat(mu, [size(imgCols,1), 1]);
    
    end
    halfSize = ceil(filterSize(1)/2);
  
    %% Low-level feature extraction.
    if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
        for filtItr = 1:filterCount
            currentFilter = double(options.filters{filtItr});

                responseImg = conv2(img, currentFilter, 'same');

            % Save response for future processing
            responseImgs(:,:,filtItr) = responseImg;
        end
    else
        % Pre-process the image blocks by whitening them.
        imgCols2 = imgCols - muArr;
        imgCols2 = imgCols2 * whMat;

        % Instead of convolving imgCols2 with the filter, we simply
        % assign a label(filter id) to each row in imgCols2 by
        % finding the cluster with minimum distance to each sample, 
        % given in every row. We also allow a soft-competition of
        % filters by removing roughly half of them, based on the
        % average distance of cluster centers from each sample.
        numberOfCols = size(imgCols2,1);
        imgCols2 = imgCols2(floor((0:((numberOfCols^2)-1)) / numberOfCols) + 1, :);
        repFilterMatrix = repmat(filterMatrix, numberOfCols, 1);
        
        % Subtract repFilterMatrix from imgCols2.
        colFiltDistances = sum((imgCols2 - repFilterMatrix).^2,2);
        distancesPerCol = reshape(colFiltDistances, numberOfCols, filterCount)';
        
        % Find average of every row, and subtract its mean from every row.
        colMeans = repmat(mean(distancesPerCol,2), 1, filterCount);
        
        % Suppress distances for every row which is more than its average
        % distance to every filter.
        distancesPerCol(distancesPerCol > colMeans) = max(max(distancesPerCol));
        
        % Find responses by reversing distances into normalized similarities.
        responses = (max(max(distancesPerCol)) - distancesPerCol) / ...
            max(max(distancesPerCol));
        
        % Assign responses to the actual image.
        realCoordIdx1 = idx1 + halfSize - 1;
        realCoordIdx2 = idx2 + halfSize - 1;
        responseImgs(realCoordIdx1, realCoordIdx2, :) = reshape(responses, [numel(idx1), numel(idx2), filterCount]);
        
        % Remove responses resulting from dead features. (Will be replaced 
        % by a better method in the future)
        deadFeatures = find(std(filterMatrix, 0, 2) < deadFeatureStd );
        for deadFeature = deadFeatures
            responseImgs(:,:,deadFeature) = 0;
        end
    end
    
    %% In Gabor-based features, we apply a minimum response threshold over response image.
    if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
       gaborFilterThr = options.gaborFilterThr * max(max(max(responseImgs)));
       responseImgs(responseImgs<max(gaborFilterThr, options.absGaborFilterThr)) = 0;
    elseif strcmp(options.filterType, 'auto')
       filterThr = options.autoFilterThr * max(max(max(responseImgs)));
       responseImgs(responseImgs<filterThr) = 0;
    end
    
   %% Write smooth object boundaries to an image based on responseImgs.
    smoothedImg = max(responseImgs,[],3);
    smoothedImg = (smoothedImg - min(min(smoothedImg))) / (max(max(smoothedImg)) - min(min(smoothedImg)));
   
   %% Inhibit weak responses in vicinity of powerful peaks.
   if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
        inhibitionHalfSize = options.gabor.inhibitionRadius;
        filterSize = options.gaborFilterSize;
   else
        filterSize = options.autoFilterSize; 
        inhibitionHalfSize = options.auto.inhibitionRadius;
   end
   halfSize=floor(filterSize/2);
   responseImgs([1:inhibitionHalfSize, (end-inhibitionHalfSize):end],:) = 0;
   responseImgs(:,[1:inhibitionHalfSize, (end-inhibitionHalfSize):end]) = 0;
   
   %% Here, we will run a loop till we clear all weak responses.
   peaks = find(responseImgs);
   peakCount = numel(peaks);
   [~, orderedPeakIdx] = sort(responseImgs(peaks), 'descend');
   orderedPeaks = peaks(orderedPeakIdx);
   validPeaks = ones(size(orderedPeaks))>0;
   [xInd, yInd, ~] = ind2sub(size(responseImgs), orderedPeaks);
   for peakItr = 1:(peakCount-1)
       if validPeaks(peakItr)
           nextPeakItr = peakItr+1;
           nearbyPeakIdx = ~(xInd(nextPeakItr:end) >= (xInd(peakItr) - inhibitionHalfSize) & xInd(nextPeakItr:end) <= (xInd(peakItr) + inhibitionHalfSize) & ...
                yInd(nextPeakItr:end) >= (yInd(peakItr) - inhibitionHalfSize) & yInd(nextPeakItr:end) <= (yInd(peakItr) + inhibitionHalfSize));
           validPeaks(nextPeakItr:end) = nearbyPeakIdx & validPeaks(nextPeakItr:end);
       end
   end
   responseImgs(orderedPeaks(~validPeaks)) = 0;
   
   % Write the responses in the final image.
   responseImgs = double(responseImgs>0);
   for filtItr = 1:filterCount
      responseImgs(:,:,filtItr) = responseImgs(:,:,filtItr) .* filtItr;
   end
   responseImg = sum(responseImgs,3);
   responseImg([1:halfSize, (end-halfSize):end],:) = 0;
   responseImg(:,[1:halfSize, (end-halfSize):end]) = 0;
   
   %% Eliminate nodes outside GT mask. If gt is not used, this does not have effect.
   responseImg(~gtMask) = 0;
   
   %% Out of this response image, we will create the nodes and output them.
   finalNodeIdx = find(responseImg);
   nodes = cell(numel(finalNodeIdx), 2);
   for nodeItr = 1:numel(finalNodeIdx)
       [centerX, centerY] = ind2sub(size(responseImg), finalNodeIdx(nodeItr));
       nodes(nodeItr,:) = {responseImg(finalNodeIdx(nodeItr)), round([centerX, centerY])};
   end
end