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
function [ nodes, activationImg, nodeActivations, smoothActivationImg ] = getNodes( img, gtFileName, options )
    %% Step 1: Get grayscaled image and assign method parameters.
    if strcmp(options.filterType, 'gabor')
        stride = options.gabor.stride;
        if size(img,3)>1
            img = rgb2gray(img(:,:,1:3));
        end
        
        %% Apply denoising to get better responses.
        for bandItr = 1:size(img,3);
    %        myfilter = fspecial('gaussian',[3 3], 2);
    %        img(:,:,bandItr) = imfilter(img(:,:,bandItr), myfilter, 'replicate', 'same', 'conv');
            img(:,:,bandItr)=medfilt2(img(:,:,bandItr), [3,3]);
        end
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
    borderSize = ceil(size(filterSize,1) / stride) + 1;
    
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
        for filtItr = 1:filterCount
            currentFilter = double(options.filters{filtItr});
            responseImg = conv2(img, currentFilter, 'same');

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
        if numberOfCols>1000 % Depends on feature dimension, but 1000 should 
                             % be reasonable for features smaller than 15x15x3)
            numberOfSets = ceil(numberOfCols/1000);
            smallRepFilterMatrix = repmat(filterMatrix, 1000, 1);
            colFiltDistances = zeros(numberOfCols * filterCount,1);
            for setItr = 1:numberOfSets
                colsInSet = min(1000, (numberOfCols-((setItr-1)*1000)));
                setCols = imgCols2((((setItr-1) * 1000) + 1):min(setItr*1000, numberOfCols),:);
                setCols = setCols(floor((0:((colsInSet * filterCount)-1)) / filterCount) + 1, :);
                if colsInSet ~= 1000
                    smallRepFilterMatrix = smallRepFilterMatrix(1:(colsInSet * filterCount),:);
                end

                % Subtract repFilterMatrix from imgCols2.
                totalAssgnCount = colsInSet * filterCount;
                assgnStartIdx = (1000 * filterCount * (setItr-1)) + 1;
                colFiltDistances(assgnStartIdx:...
                    (assgnStartIdx+totalAssgnCount-1)) = ...
                    sqrt(sum((setCols - smallRepFilterMatrix).^2,2));
            end
            clear smallRepFilterMatrix imgCols2 setCols;
        else
            % No need to divide, just find the distances.
            imgCols2 = imgCols2(floor((0:((numberOfCols * filterCount)-1)) / filterCount) + 1, :);
            repFilterMatrix = repmat(filterMatrix, numberOfCols, 1);

            % Subtract repFilterMatrix from imgCols2.
            colFiltDistances = sqrt(sum((imgCols2 - repFilterMatrix).^2,2));
            clear repFilterMatrix imgCols2;
        end
        clear imgCols;
        % Reshape distances into a NxD array, where N is number of columns (blocks),
        % and D is number of filters.
        distancesPerCol = reshape(colFiltDistances, filterCount, numberOfCols).';
        
        % Find average of every row.
        colMeans = repmat(mean(distancesPerCol,2), 1, filterCount);
        
        % Suppress distances for every row which is more than its average
        % distance to every filter.
        meanAssgnIdx = distancesPerCol > colMeans;
        distancesPerCol(meanAssgnIdx) = colMeans(meanAssgnIdx);
        
        % Find responses by reversing distances into normalized similarities.
%         responses = (colMeans - distancesPerCol) ./ ...
%             colMeans;
        responses = (colMeans - distancesPerCol);
        
        % Assign responses to the actual image.
        realCoordIdx1 = idx1 + halfSize - 1;
        realCoordIdx2 = idx2 + halfSize - 1;
        [p, q] = meshgrid(realCoordIdx2, realCoordIdx1);
        realCoords = [q(:), p(:)];
        responseImgs = reshape(responses, [numel(idx1), numel(idx2), filterCount]);
        
        % Subsample gt mask.
        gtMask = imresize(gtMask((halfSize+1):(end-halfSize), (halfSize+1):(end-halfSize)), newImgSize);
    end
    
    %% We apply a minimum response threshold over response image.
    if strcmp(options.filterType, 'gabor')
       filterThr = options.gaborFilterThr * max(max(max(responseImgs)));
       responseImgs(responseImgs<max(filterThr, options.absGaborFilterThr)) = 0;
    else
       filterThr = options.autoFilterThr * max(max(max(responseImgs)));
       responseImgs(responseImgs<filterThr) = 0;
    end
    
%   %% Write smooth object boundaries to an image based on responseImgs.
    smoothActivationImg = max(responseImgs,[],3);
    smoothActivationImg = (smoothActivationImg - min(min(smoothActivationImg))) / (max(max(smoothActivationImg)) - min(min(smoothActivationImg)));
   
   %% Inhibit weak responses in vicinity of powerful peaks.
    if strcmp(options.filterType, 'gabor')
        inhibitionHalfSize = options.gabor.inhibitionRadius;
    else
        inhibitionHalfSize = options.auto.inhibitionRadius;
    end
    responseImgs([1:inhibitionHalfSize, (end-inhibitionHalfSize):end],:) = 0;
    responseImgs(:,[1:inhibitionHalfSize, (end-inhibitionHalfSize):end]) = 0;

    % Each response will clear other weak responses at the very same pixel.
    % Use this feature to get rid of most peaks.
    [activationImg, nodeIdImg] = max(responseImgs, [], 3);
    clear responseImgs;
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
    responseImg([1:borderSize, (end-borderSize):end],:) = 0;
    responseImg(:,[1:borderSize, (end-borderSize):end]) = 0;
    activationImg([1:borderSize, (end-borderSize):end],:) = 0;
    activationImg(:,[1:borderSize, (end-borderSize):end]) = 0;
    activationImg = activationImg / max(max(activationImg));

    %% Eliminate nodes outside GT mask. If gt is not used, this does not have effect.
    responseImg(~gtMask) = 0;

    %% Out of this response image, we will create the nodes and output them.
    responseImg(ismember(responseImg, deadFeatures)) = 0;
    activationImg(responseImg == 0) = 0;
    finalNodeIdx = find(responseImg);
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
       [centerX, centerY] = ind2sub(size(responseImg), finalNodeIdx(nodeItr));
       nodes(nodeItr,:) = {responseImg(finalNodeIdx(nodeItr)), round([centerX, centerY]), round(realCoords(nodeItr,:))};
    end
    nodes = nodes(cellfun(@(x) ~isempty(x), nodes(:,1)),:);
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