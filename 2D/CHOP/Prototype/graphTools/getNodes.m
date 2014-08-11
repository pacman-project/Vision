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
        whMat = options.auto.whMat;
        invMat = options.auto.invMat;
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
    tempImgs = zeros(size(img,1), size(img,2), size(img,3));
  %  filterBand = [];
    imgCols = zeros((size(img,1)-filterSize(1)+1) * (size(img,2)-filterSize(1)+1), prod(filterSize));
    startIdx = 1;
    iterator = prod(filterBandSize)-1;
    for bandItr = 1:size(img,3)
%         filterBand = currentFilter(:,:,bandItr);
        imgCols(:,startIdx:(startIdx+iterator)) = im2col(img(:,:,bandItr), filterBandSize)';
        startIdx = startIdx + iterator + 1;
    end
    muArr = repmat(mu, [size(imgCols,1), 1]);
    halfSize = ceil(filterSize(1)/2);
  
    for filtItr = 1:filterCount
        currentFilter = double(options.filters{filtItr});
  
        if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
            responseImg = conv2(img, currentFilter, 'same');
        else
%            responseImg = mean(convn(img, currentFilter, 'same'), 3);
            imgCols2 = imgCols - muArr;
            imgCols2 = imgCols2 * whMat;
            imgCols2 = imgCols2 * currentFilter(:);
            responses = (imgCols2-min(imgCols2)) / (max(imgCols2) - min(imgCols2));
            responseImg = zeros(size(img,1), size(img,2));
            responseImg(halfSize:(end-halfSize+1), halfSize:(end-halfSize+1)) = reshape(responses, ...
                [size(responseImg,1)-filterSize(1)+1, size(responseImg,2)-filterSize(1)+1]);
        end
            
        % Save response for future processing
        responseImgs(:,:,filtItr) = responseImg;
    end
    
%     function x = myConv(x)
%         x = x - mu;
%         x = x * whMat;
%         x = sum(x*currentFilter);
%     end
    
    %% In Gabor-based features, we apply a minimum response threshold over response image.
    if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
       gaborFilterThr = options.gaborFilterThr * max(max(max(responseImgs)));
    %    gaborFilterThr = 150;
       responseImgs(responseImgs<max(gaborFilterThr, options.absGaborFilterThr)) = 0;
    %     responseImgs(responseImgs<gaborFilterThr) = 0;
    elseif strcmp(options.filterType, 'auto')
       filterThr = options.autoFilterThr * max(max(max(responseImgs)));
       responseImgs(responseImgs<filterThr) = 0;
    end
    
   %% Write smooth object boundaries to an image based on responseImgs.
%   smoothedImg = getSmoothedImage(responseImgs, options.filters);
    smoothedImg = mean(responseImgs,3);
    smoothedImg = (smoothedImg - min(min(smoothedImg))) / (max(max(smoothedImg)) - min(min(smoothedImg)));
   
   %% Inhibit weak responses in vicinity of powerful peaks.
   inhibitionHalfSize = options.gabor.inhibitionRadius;
   if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
      filterSize = options.gaborFilterSize;
   else
      filterSize = options.autoFilterSize; 
   end
   halfSize=floor(filterSize/2);
   responseImgs([1:inhibitionHalfSize, (end-inhibitionHalfSize):end],:) = 0;
   responseImgs(:,[1:inhibitionHalfSize, (end-inhibitionHalfSize):end]) = 0;
   
   %% Here, we will run a loop till we clear all weak responses.
   peaks = find(responseImgs);
%   prevPeakCount = inf;
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
%   peaks = find(validPeaks);
%   peakCount = numel(peaks);
%    
%    while prevPeakCount ~= peakCount
%        [xInd, yInd, ~] = ind2sub(size(responseImgs), peaks);
%        for peakItr = 1:size(peaks,1)
%           % If this peak has not yet been eliminated, go check nearby peaks
%           if responseImgs(peaks(peakItr)) > 0 
%               nearbyPeakIdx = xInd >= (xInd(peakItr) - inhibitionHalfSize) & xInd <= (xInd(peakItr) + inhibitionHalfSize) & ...
%                   yInd >= (yInd(peakItr) - inhibitionHalfSize) & yInd <= (yInd(peakItr) + inhibitionHalfSize);
% 
%               maxPeak = max(responseImgs(peaks(nearbyPeakIdx)));
%               if responseImgs(peaks(peakItr)) > (maxPeak-0.0001)
%                  nearbyPeakIdx(peakItr) = 0;
%                  responseImgs(peaks(nearbyPeakIdx)) = 0; 
%               else
%        %          responseImgs(peaks(peakItr)) = 0;
%               end
%           end
%        end
%        prevPeakCount = peakCount;
%        peaks = find(responseImgs);
%        peakCount = numel(peaks);
%    end
   
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