%> Name: GetGaborResponses( img, options )
%>
%> Description: Given an input image and gabor filters, this function
%> processes 'img' in order to obtain individual response images for the
%> convolution of each filter. 
%>
%> @param img Input image
%> @param options Program options.
%> 
%> @retval responseImgs
%>
%> Author: rusen
%>
%> Updates
%> Ver 1.0 on 25.10.2015
function [ nodes, activationImg, nodeActivations, smoothActivationImg ] = GetGaborResponses( img, options )
    %% Step 1: Get grayscaled image and assign method parameters.
   stride = options.gabor.stride;
   if size(img,3)>1
       img = rgb2gray(img(:,:,1:3));
   end

   %% Apply denoising to get better responses.
   for bandItr = 1:size(img,3);
       img(:,:,bandItr)=medfilt2(img(:,:,bandItr), [3,3]);
   end
    filterSize = size(options.filters{1});
    img = double(img);
    filterCount = numel(options.filters);
    borderSize = ceil(size(filterSize,1) / stride) + 1;
    
    %% In case of auto-learned features, get each patch (using stride) as a separate column of data.
    newImgSize = floor([(size(img,1)-1)/stride, (size(img,2)-1)/stride]) + 1;
    responseImgs = zeros(newImgSize(1), newImgSize(2), filterCount);
    
    %% Low-level feature extraction.
   realCoordIdx = zeros(filterCount, prod(newImgSize));
   for filtItr = 1:filterCount
       currentFilter = double(options.filters{filtItr});
       responseImg = conv2(img, rot90(currentFilter,2), 'same');

       % Simulate stride here, and subsample the image.
       [responseImg, idx] = MaxPooling(responseImg, [stride, stride]);
       realCoordIdx(filtItr,:) = idx';
       % Save response for future processing
       responseImgs(:,:,filtItr) = responseImg;
   end
    
    %% We apply a minimum response threshold over response image.
   filterThr = options.gaborFilterThr * max(max(max(responseImgs)));
   responseImgs(responseImgs<max(filterThr, options.absGaborFilterThr)) = 0;
    
   %% Inhibit weak responses in vicinity of powerful peaks.
    inhibitionHalfSize = options.gabor.inhibitionRadius;
    responseImgs([1:inhibitionHalfSize, ((end-inhibitionHalfSize)+1):end],:, :) = 0;
    responseImgs(:,[1:inhibitionHalfSize, ((end-inhibitionHalfSize)+1):end], :) = 0;
 
    % Each response will clear other weak responses at the very same pixel.
    % Use this feature to get rid of most peaks.
    [activationImg, nodeIdImg] = max(responseImgs, [], 3);
    smoothActivationImg = activationImg;
%    smoothActivationImg = smoothActivationImg>0;
    smoothActivationImg = smoothActivationImg/max(max(max(smoothActivationImg)));
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
    idx = sub2ind(size(realCoordIdx), responseImg(finalNodeIdx), finalNodeIdx);
    realCoordLin = realCoordIdx(idx);
    [realCoordX, realCoordY] = ind2sub(size(img), realCoordLin);
    realCoords = [realCoordX, realCoordY];
    nodes = cell(numel(finalNodeIdx), 3);
    nodeActivations = single(activationImg(finalNodeIdx));
    for nodeItr = 1:numel(finalNodeIdx)
       [centerX, centerY] = ind2sub(size(responseImg), finalNodeIdx(nodeItr));
       nodes(nodeItr,:) = {responseImg(finalNodeIdx(nodeItr)), round([centerX, centerY]), round(realCoords(nodeItr,:))};
    end
    nodes = nodes(cellfun(@(x) ~isempty(x), nodes(:,1)),:);
end