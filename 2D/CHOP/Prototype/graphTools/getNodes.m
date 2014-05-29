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
    %% Step 1: Get grayscaled and binarized image.
    rgbImg = img;
    if size(img,3)>1
        img = rgb2gray(img(:,:,1:3));
    end
    doubleImg = uint8(img);
    
    % Either binarize the image using Otsu, or extract edges using Canny.
    %edgeImg = double(edge(doubleImg(:,:), 'canny'));
%    edgeImg = im2bw(doubleImg, graythresh(doubleImg));
    edgeImg = img;
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
        gtMask = ones(size(doubleImg)) > 0;
    end
    
%     w     = 2;       % bilateral filter half-width
%     sigma = [3 0.1]; % bilateral filter standard deviations
%     edgeImg = bfilter2(double(edgeImg)/double(max(max(edgeImg))),w,sigma);
%     edgeImg = uint8(edgeImg * 255);

    % Apply filtering to get better responses.
     myfilter = fspecial('gaussian',[3 3], 2);
     edgeImg = imfilter(edgeImg, myfilter, 'replicate', 'same', 'conv');
     edgeImg=medfilt2(edgeImg, [3,3]);
%     limits = stretchlim(edgeImg);
%     edgeImg = imadjust(edgeImg, limits, []);
    
 %   limits = stretchlim(edgeImg, 0.05);
 %   edgeImg = imadjust(edgeImg, limits, []);
    
%     H = padarray(1,[2 2]) - fspecial('gaussian',[5 5], 1); % create unsharp mask
%     myfilteredimage2 = imfilter(img, H, 'replicate', 'same', 'conv');
    
    %% Get response by applying each filter to the image.
    responseImgs = zeros(size(edgeImg,1), size(edgeImg,2), filterCount);
    for filtItr = 1:filterCount
%         responseImg = zeros(size(edgeImg));
%         for bandItr = 1:size(rgbImg,3)
%             responseImg = max(responseImg, conv2(double(rgbImg(:,:,bandItr)), options.filters{filtItr}, 'same'));
%         end
       responseImg = conv2(double(edgeImg), options.filters{filtItr}, 'same');
       % Save response for future processing
       responseImgs(:,:,filtItr) = responseImg;
    end
    
   %% In Gabor-based features, we apply a minimum response threshold over response image.
   if strcmp(options.filterType, 'gabor') || strcmp(options.filterType, 'lhop')
       gaborFilterThr = options.gaborFilterThr * max(max(max(responseImgs)));
   %    gaborFilterThr = 150;
       responseImgs(responseImgs<max(gaborFilterThr, options.absGaborFilterThr)) = 0;
   %     responseImgs(responseImgs<gaborFilterThr) = 0;
   end
    
   %% Write smooth object boundaries to an image based on responseImgs.
%   smoothedImg = getSmoothedImage(responseImgs, options.filters);
    smoothedImg = mean(responseImgs,3);
    smoothedImg = (smoothedImg - min(min(smoothedImg))) / (max(max(smoothedImg)) - min(min(smoothedImg)));
   
   %% Inhibit weak responses in vicinity of powerful peaks.
   inhibitionHalfSize = options.gabor.inhibitionRadius;
   halfSize=floor(options.gaborFilterSize/2);
   responseImgs([1:inhibitionHalfSize, (end-inhibitionHalfSize):end],:) = 0;
   responseImgs(:,[1:inhibitionHalfSize, (end-inhibitionHalfSize):end]) = 0;
   
   %% Here, we will run a loop till we clear all weak responses.
   peaks = find(responseImgs);
   prevPeakCount = inf;
   peakCount = numel(peaks);
   while prevPeakCount ~= peakCount
       [xInd, yInd, ~] = ind2sub(size(responseImgs), peaks);
       for peakItr = 1:size(peaks,1)
          % If this peak has not yet been eliminated, go check nearby peaks
          if responseImgs(peaks(peakItr)) > 0 
              nearbyPeakIdx = xInd >= (xInd(peakItr) - inhibitionHalfSize) & xInd <= (xInd(peakItr) + inhibitionHalfSize) & ...
                  yInd >= (yInd(peakItr) - inhibitionHalfSize) & yInd <= (yInd(peakItr) + inhibitionHalfSize);

              maxPeak = max(responseImgs(peaks(nearbyPeakIdx)));
              if responseImgs(peaks(peakItr)) > (maxPeak-0.0001)
                 nearbyPeakIdx(peakItr) = 0;
                 responseImgs(peaks(nearbyPeakIdx)) = 0; 
              else
       %          responseImgs(peaks(peakItr)) = 0;
              end
          end
       end
       prevPeakCount = peakCount;
       peaks = find(responseImgs);
       peakCount = numel(peaks);
   end
   
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

%> Name: getSmoothedImage
%>
%> Description: Writes each filter response to the output image to form a
%> nearly correct reconstruction based on current responses.
%>
%> @param responseImgs Response images with the same number as filters.
%> @param filters The filters applied.
%>
%> @retval smoothedImg The output response.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 12.01.2014
function [smoothedImg] = getSmoothedImage(responseImgs, filters)
    smoothedImg = zeros(size(responseImgs,1), size(responseImgs,2));
    halfSize = floor(size(filters{1},1)/2);
    responseImgs([1:(halfSize+1), (end-halfSize):end], :, :) = 0;
    responseImgs(:, [1:(halfSize+1), (end-halfSize):end], :) = 0;
    responseImgs = responseImgs/max(max(max(responseImgs)));
    for filterItr = 1:size(filters,1)
        [xInd, yInd, val] = find(responseImgs(:,:,filterItr));
        for peakItr = 1:size(val,1)
            smoothedImg((xInd(peakItr)-halfSize):(xInd(peakItr)+halfSize), ...
                (yInd(peakItr)-halfSize):(yInd(peakItr)+halfSize)) = val(peakItr) * filters{filterItr};
        end
        smoothedImg = smoothedImg + conv2(responseImgs(:,:,filterItr), filters{filterItr}, 'same');
    end
    smoothedImg = smoothedImg/max(max(smoothedImg));
end

