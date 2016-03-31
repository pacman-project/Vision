%> Name: getInitialNodes
%>
%> Description: When given an image, this function is used to get initial
%> simple nodes by applying low-level filters over edge mask of the image. Each
%> response peak is considered as a simple feature instance.
%>
%> @param img Input image
%> @param gtFileName An empty param either means the gt for that image is not
%> given, or it is of invalid size and cannot be used. If not empty, we
%> use it in elimination of nodes (given the gt use is enabled in 
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
function [ nodes, activationImg, nodeActivations, smoothActivationImg, trueResponseImgs, edgeImg ] = getNodes( img, gtFileName, options )
    % Step 0: Obtain background
    if size(img,3)>1
         grayImg = rgb2gray(img(:,:,1:3));
    else
         grayImg = img;
    end
    backgroundMask = grayImg == 0;
    backgroundMask = imdilate(backgroundMask, strel('disk',3));
    edgeMask = edge(img, 'canny', [0.05 0.1], 1);
    outEdgeMask = backgroundMask & edgeMask;
    inEdgeMask = ~backgroundMask & edgeMask;
    edgeImg = double(outEdgeMask) + double(inEdgeMask) * 0.5;
%     edgeMask = imdilate(edgeMask,strel('disk', 2));
%     edgeMask(backgroundMask) = 0;
%     
    % In addition, we filter out weak responses for outside borders of the
    % object.

    %% Step 1: Get grayscaled image and assign method parameters.
    if strcmp(options.filterType, 'gabor')
        stride = options.gabor.stride;
        if size(img,3)>1
            img = rgb2gray(img(:,:,1:3));
        end
        
        %% Apply denoising to get better responses.
%         for bandItr = 1:size(img,3);
%     %        myfilter = fspecial('gaussian',[3 3], 2);
%     %        img(:,:,bandItr) = imfilter(img(:,:,bandItr), myfilter, 'replicate', 'same', 'conv');
%             img(:,:,bandItr)=medfilt2(img(:,:,bandItr), [3,3]);
%         end
        deadFeatures = [];
    elseif strcmp(options.filterType, 'auto')
        stride = options.auto.stride;
        whMat = options.auto.whMat;
        mu = options.auto.mu;
        deadFeatures = options.auto.deadFeatures;
        filterMatrix = options.filterMatrix;
    else
        nodes = [];
        display('Feature type not implemented (in getNodes.m).');
        return;
    end
    filterSize = size(options.filters{1});
    filterBandSize = filterSize(1:2);
    img = double(img);
    filterCount = numel(options.filters);
    
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
    
    %% In case of auto-learned features, get each patch (using stride) as a separate column of data.
    newImgSize = floor([(size(img,1)-1)/stride, (size(img,2)-1)/stride]) + 1;
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
        newImgSize = [numel(idx1), numel(idx2)];
        
        for bandItr = 1:size(img,3)
            tempCols = im2col(img(:,:,bandItr), filterBandSize)';
            imgCols(:,startIdx:(startIdx+iterator)) = tempCols(validCols,:);
            startIdx = startIdx + iterator + 1;
        end
        clear tempCols;
        muArr = repmat(mu, [size(imgCols,1), 1]);
    end
    responseImgs = zeros(newImgSize(1), newImgSize(2), filterCount);
    halfSize = ceil(filterSize(1)/2);
    
    %% Low-level feature extraction.
    if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
        realCoordIdx = zeros(filterCount, prod(newImgSize));
        
        % We're using steerable filters here. Mimicking gabor filters.
%        [res, theta, nms, rotations] = steerableDetector(img, 5, 0.3, filterCount * 2);
          [~, ~, ~, rotations] = steerableDetector(img, 5, 1, filterCount * 2);
 %       filterOrder = [1,4,3,2,5,8,7,6,9,12,11,10,13,16,15,14];
%        rotations = rotations(:,:,filterOrder);

        trueResponseImgs = zeros(size(img,1), size(img,2), filterCount);
        outMaxVal = 1;
        inMaxVal = 1;
        for filtItr = 1:filterCount
%              filtId = rem(filterCount+2-filtItr, filterCount);
%              if filtId == 0
%                   filtId = filtId + filterCount;
%              end
 %           currentFilter = double(options.filters{filtItr});
 %           responseImg = conv2(img, rot90(currentFilter,2), 'same');
            responseImg = max(rotations(:,:,filtItr), rotations(:,:,filtItr+filterCount));
%            
%            figure, imshow(responseImg);
            % Process response image so it only exists under edges.
            outMaxVal = max(outMaxVal, max(max(responseImg(outEdgeMask))));
            inMaxVal = max(inMaxVal, max(max(responseImg(inEdgeMask))));
            trueResponseImgs(:,:,filtItr) = responseImg;
        end
        if isempty(inMaxVal)
             inMaxVal = 1;
        end
        if isempty(outMaxVal)
             outMaxVal = 1;
        end
       filterThr = options.gaborFilterThr * outMaxVal;
       inFilterThr = options.innerGaborFilterThr * outMaxVal;
       
       [smoothActivationImg, ~] = max(trueResponseImgs, [], 3);
       smoothActivationImg = smoothActivationImg / max(max(smoothActivationImg));
        
        % Obtain a smooth activation image.
        for filtItr = 1:filterCount
            responseImg = squeeze(trueResponseImgs(:,:,filtItr));
            
            % Apply filter thresholds.
            if ~isempty(filterThr)
                responseImg(outEdgeMask) = responseImg(outEdgeMask) / outMaxVal;
            end
            if ~isempty(inFilterThr)
                responseImg(inEdgeMask) = responseImg(inEdgeMask) / inMaxVal;
            end
            responseImg(responseImg>1) = 1;
            trueResponseImgs(:,:,filtItr) = responseImg;
       end
        
       % Eliminate weak responses and perform pooling.
       for filtItr = 1:filterCount
            responseImg = squeeze(trueResponseImgs(:,:,filtItr));
            responseImg(~edgeMask) = 0;
            
            % Apply filter thresholds.
            if ~isempty(filterThr)
                responseImg(outEdgeMask) = double(responseImg(outEdgeMask) * outMaxVal >= filterThr) .* responseImg(outEdgeMask);
            end
            if ~isempty(inFilterThr)
                responseImg(inEdgeMask) = double(responseImg(inEdgeMask) * inMaxVal >= inFilterThr) .* responseImg(inEdgeMask);
            end
            % Simulate stride here, and subsample the image.
            [responseImg, idx] = MaxPooling(responseImg, [stride, stride]);
            realCoordIdx(filtItr,:) = idx';
            % Save response for future processing
            responseImgs(:,:,filtItr) = responseImg;
       end
        
        % Subsample gt mask.
        gtMask = imresize(gtMask, newImgSize);
    else
        % Pre-process the image blocks by whitening them.
        if strcmp(options.autoNormalize, 'whiten')
             imgCols2 = imgCols - muArr;
             imgCols2 = imgCols2 * whMat;
        else
             imgCols2 = normalizeData(imgCols);
        end

        % Instead of convolving imgCols2 with the filter, we simply
        % assign a label(filter id) to each row in imgCols2 by
        % finding the cluster with minimum distance to each sample, 
        % given in every row. We also allow a soft-competition of
        % filters by removing roughly half of them, based on the
        % average distance of cluster centers from each sample.
        numberOfCols = size(imgCols2,1);
        % If number of cols is more than 1000, we divide and conquer
        % imgCols2 array into fixed-size bins. Otherwise, it's hell of a
        % memory problem.
%        if numberOfCols>1000 % Depends on feature dimension, but 1000 should 
%                              % be reasonable for features smaller than 15x15x3)
%             numberOfSets = ceil(numberOfCols/1000);
%             smallRepFilterMatrix = repmat(filterMatrix, 1000, 1);
%             colFiltDistances = zeros(numberOfCols * filterCount,1);
%             for setItr = 1:numberOfSets
%                 colsInSet = min(1000, (numberOfCols-((setItr-1)*1000)));
%                 setCols = imgCols2((((setItr-1) * 1000) + 1):min(setItr*1000, numberOfCols),:);
%                 setCols = setCols(floor((0:((colsInSet * filterCount)-1)) / filterCount) + 1, :);
%                 if colsInSet ~= 1000
%                     smallRepFilterMatrix = smallRepFilterMatrix(1:(colsInSet * filterCount),:);
%                 end
% 
%                 % Subtract repFilterMatrix from imgCols2.
%                 totalAssgnCount = colsInSet * filterCount;
%                 assgnStartIdx = (1000 * filterCount * (setItr-1)) + 1;
%                 colFiltDistances(assgnStartIdx:...
%                     (assgnStartIdx+totalAssgnCount-1)) = ...
%                     sqrt(sum((setCols - smallRepFilterMatrix).^2,2));
%             end
%             clear smallRepFilterMatrix imgCols2 setCols;
%        else
        % No need to divide, just find the distances.
        responses = zeros(numberOfCols, filterCount);
        for filterItr = 1:filterCount
            responses(:,filterItr) = sqrt(sum((imgCols2 - repmat(filterMatrix(filterItr,:), numberOfCols,1)).^2, 2));
%            responses(:,filterItr) = responses(:,filterItr) / max(responses(:,filterItr));
        end
        responses = max(max(responses)) - responses;
        responses = responses / max(max(max(responses)));
        clear imgCols imgCols2;
        
        % Assign responses to the actual image.
        realCoordIdx1 = idx1 + halfSize - 1;
        realCoordIdx2 = idx2 + halfSize - 1;
        [p, q] = meshgrid(realCoordIdx2, realCoordIdx1);
        realCoords = [q(:), p(:)];
        responseImgs = reshape(responses, [numel(idx1), numel(idx2), filterCount]);
        clear responses realCoordIdx2;
        
        % Subsample gt mask.
        gtMask = imresize(gtMask((halfSize+1):(end-halfSize), (halfSize+1):(end-halfSize)), newImgSize);
    end
    
    %% We apply a minimum response threshold over response image.
    % We treat background-foreground edges differently.
%     if strcmp(options.filterType, 'gabor')
%        filterThr = max(options.absGaborFilterThr, options.gaborFilterThr * max(max(max(responseImgs))));
%        inFilterThr = options.innerGaborFilterThr * max(max(max(responseImgs)));
%     else
%        filterThr = max(options.autoFilterThr, options.autoFilterThr * max(max(max(responseImgs))));
%     end
%     if strcmp(options.filterType, 'gabor')
%          backgroundMask = imresize(backgroundMask, [size(responseImgs,1), size(responseImgs,2)], 'bilinear');
%          filterMask = ~backgroundMask & imdilate(backgroundMask, strel('disk', 2));
%          innerMask = ~filterMask & ~backgroundMask;
%          maxOutVal = 1;
%          maxInVal = 1;
%          
%          % Find maximum in/out response values.
%          for bandItr = 1:size(responseImgs,3)
%               % Obtain responses.
%               tempImg = squeeze(responseImgs(:,:,bandItr));
%               
%               % Process outer edges.
%               maxOutVal = max(maxOutVal, max(tempImg(filterMask)));
%               
%               % Process inner edges.
%               maxInVal = max(maxInVal, max(tempImg(innerMask)));
%          end
%          
%          for bandItr = 1:size(responseImgs,3)
%               % Obtain responses.
%               tempImg = squeeze(responseImgs(:,:,bandItr));
%               
%               % Process outer edges.
%               vals = tempImg(filterMask);
%               vals(vals<filterThr) = 0;
%               tempImg(filterMask) = vals;
%               tempImg(filterMask) = vals / maxOutVal;
%               
%               % Process inner edges.
%               vals = tempImg(innerMask);
%               vals(vals < inFilterThr) = 0;
%               tempImg(innerMask) = vals;
%               tempImg(innerMask) = vals / maxInVal;
%               
%               % Set outside points to zero.
%               tempImg(backgroundMask) = 0;
%               
%               % Write everything back.
%               responseImgs(:,:,bandItr) = tempImg;
%          end
%     else
%         responseImgs(responseImgs<filterThr) = 0;
%     end
    
   %% Inhibit weak responses in vicinity of powerful peaks.
    if strcmp(options.filterType, 'gabor')
        inhibitionHalfSize = options.gabor.inhibitionRadius;
    else
        inhibitionHalfSize = options.auto.inhibitionRadius;
    end
    responseImgs([1:inhibitionHalfSize, ((end-inhibitionHalfSize)+1):end],:, :) = 0;
    responseImgs(:,[1:inhibitionHalfSize, ((end-inhibitionHalfSize)+1):end], :) = 0;

    % Each response will clear other weak responses at the very same pixel.
    % Use this feature to get rid of most peaks.
    [activationImg, nodeIdImg] = max(responseImgs, [], 3);
    
    % Process only edges.
%     if strcmp(options.filterType, 'gabor')
%          edgeMask = imresize(edgeMask, size(activationImg));
%          activationImg(~edgeMask) = 0;
%          activationImg(activationImg > 0) = 1;
%          nodeIdImg(~edgeMask) = 0;
%     end
%    smoothActivationImg = activationImg;
    peaks = find(activationImg);
    
    %% Here, we will run a loop till we clear all weak responses.
    peakCount = numel(peaks);
    [~, orderedPeakIdx] = sort(activationImg(peaks), 'descend');
    orderedPeaks = peaks(orderedPeakIdx);
    validPeaks = ones(size(orderedPeaks))>0;
    
    % If inhibition is needed, run it.
    if inhibitionHalfSize > 0
         [xInd, yInd, ~] = ind2sub(size(activationImg), orderedPeaks);
         for peakItr = 1:(peakCount-1)
            if validPeaks(peakItr)
                nextPeakItr = peakItr+1;
                nearbyPeakIdx = ~(xInd(nextPeakItr:end) >= (xInd(peakItr) - inhibitionHalfSize) & xInd(nextPeakItr:end) <= (xInd(peakItr) + inhibitionHalfSize) & ...
                     yInd(nextPeakItr:end) >= (yInd(peakItr) - inhibitionHalfSize) & yInd(nextPeakItr:end) <= (yInd(peakItr) + inhibitionHalfSize));
                validPeaks(nextPeakItr:end) = nearbyPeakIdx & validPeaks(nextPeakItr:end);
            end
         end
         activationImg(orderedPeaks(~validPeaks)) = 0;
    end
    
    % Write the responses in the final image.
    responseImg = zeros(size(activationImg));
    responseImg(orderedPeaks(validPeaks)) = nodeIdImg(orderedPeaks(validPeaks));
    borderSize = ceil(halfSize / stride);
    responseImg([1:borderSize, (end+1-borderSize):end],:) = 0;
    responseImg(:,[1:borderSize, (end+1-borderSize):end]) = 0;
    activationImg([1:borderSize, (end+1-borderSize):end],:) = 0;
    activationImg(:,[1:borderSize, (end+1-borderSize):end]) = 0;
    activationImg = activationImg / max(max(activationImg));

    %% Eliminate nodes outside GT mask. If gt is not used, this does not have effect.
    responseImg(~gtMask) = 0;
    responseImg(ismember(responseImg, deadFeatures)) = 0;
    activationImg(responseImg == 0) = 0;
    
    %% We perform pooling with a stride defined by local inhibition.
%     if inhibitionHalfSize > 0
%          [vals, valsIdx] = MaxPooling(responseImg, [inhibitionHalfSize+1, inhibitionHalfSize+1]);
%          sizeVals =size(vals);
%          finalNodeIdx = sort(valsIdx(vals>0));
%     else
         finalNodeIdx = find(responseImg);
%    end
    
    %% Out of this response image, we will create the nodes and output them.
    if strcmp(options.filterType, 'gabor')
         idx = sub2ind(size(realCoordIdx), responseImg(finalNodeIdx), finalNodeIdx);
         realCoordLin = realCoordIdx(idx);
         [realCoordX, realCoordY] = ind2sub(size(img), realCoordLin);
         realCoords = [realCoordX, realCoordY];
    else
         realCoords = realCoords(finalNodeIdx,:); 
    end
    nodes = cell(numel(finalNodeIdx), 3);
    nodeActivations = single(activationImg(finalNodeIdx));
    for nodeItr = 1:numel(finalNodeIdx)
%        if inhibitionHalfSize > 0
%            centerIdx = find(valsIdx == finalNodeIdx(nodeItr));
%            [centerX, centerY] = ind2sub(sizeVals, centerIdx);
%        else
           [centerX, centerY] = ind2sub(size(responseImg), finalNodeIdx(nodeItr));
%       end
       nodes(nodeItr,:) = {responseImg(finalNodeIdx(nodeItr)), round([centerX, centerY]), round(realCoords(nodeItr,:))};
    end
    nodes = nodes(cellfun(@(x) ~isempty(x), nodes(:,1)),:);
     
    % Prepare output image.
    activationImg = zeros(size(img,1), size(img,2));
    idx = sub2ind(size(activationImg), realCoords(:,1), realCoords(:,2));
    activationImg(idx) = nodeActivations;
end

function patches = normalizeData(patches)

     % Squash data to [0.1, 0.9] since we use sigmoid as the activation
     % function in the output layer

     % Remove DC (mean of images). 
     patches = bsxfun(@minus, patches, mean(patches));

     % Truncate to +/-3 standard deviations and scale to -1 to 1
     pstd = 3 * std(patches(:));
     patches = max(min(patches, pstd), -pstd) / pstd;

     % Rescale from [-1,1] to [0.1,0.9]
     patches = (patches + 1) * 0.4 + 0.1;
end