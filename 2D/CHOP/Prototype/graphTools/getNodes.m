%> Name: getInitialNodes
%>
%> Description: When given an image, this function is used to get initial
%> simple nodes by applying Gabor filters over edge mask of the image. Each
%> response peak is considered as a simple feature instance.
%>
%> @param img Input image
%>
%> @param filterCount (default 6) Number of rotated filters to apply on the
%> edge map.
%>
%> @param gaborFilterThr Maximum response * gaborFilterThr will be accepted 
%> as a valid response. In [0,1].
%>
%> @param gaborFilterSize Size of the filter to be applied over image.
%>
%> @param gaborAreaMinResponse Minimum response in the Gabor filter to
%> count is part of the main filter. The idea is that each filter response
%> inhibits other lower-level responses in its area, specified by this 
%> threshold applied over the filter.
%>
%> @param currentFolder The workspace folder.
%> 
%> @retval nodes The nodes to form further graphs.
%>
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 18.11.2013
%> Ver 1.1 on 03.12.2013 Response inhibition added.

function [ nodes, smoothedImg ] = getNodes( img, filterCount, gaborFilterThr, gaborAreaMinResponse, gaborFilterSize, currentFolder)
    %% Get grayscaled image and get the edges/
    if size(img,3)>1
        img = rgb2gray(img(:,:,1:3));
    end
    doubleImg = img;
    
%    edgeImg = double(edge(doubleImg(:,:), 'canny'));
    edgeImg = im2bw(doubleImg, graythresh(doubleImg));

    % Run each filter on the band and collect responses
    refGaborFilt = gabor_fn(1, 0, 1.2, 0, 0.4);
    dimSize = size(refGaborFilt);
    
    %% Pad or crop the sides of the gabor filter to set it to a fixed size of gaborFilterSize x gaborFilterSize
    newFilt = zeros(gaborFilterSize);
    if dimSize(1) >= gaborFilterSize
       lowerPad = ceil((dimSize(1)-gaborFilterSize)/2);
       upperPad = (dimSize(1) - lowerPad) - gaborFilterSize;
       refGaborFilt = refGaborFilt((lowerPad+1):(end-upperPad), :);
    end
    if dimSize(2) >= gaborFilterSize
       leftPad = ceil((dimSize(2)-gaborFilterSize)/2);
       rightPad = (dimSize(2) - leftPad) - gaborFilterSize;
       refGaborFilt = refGaborFilt((leftPad+1):(end-rightPad), :);
    end
    dimSize = size(refGaborFilt);
    lowerPad = ceil((gaborFilterSize - dimSize(1))/2);
    upperPad = (gaborFilterSize - lowerPad) - dimSize(1);
    leftPad = ceil((gaborFilterSize - dimSize(2))/2);
    rightPad = (gaborFilterSize - leftPad) - dimSize(2);
    newFilt((lowerPad+1):(end-upperPad), ((leftPad+1):(end-rightPad))) = refGaborFilt;
    refGaborFilt = newFilt;
    
    % Learn normalization constant
    normConst = sum(sum(refGaborFilt));
    
    %% Get response by rotating the reference filter and applying over image.
    responseImgs = zeros(size(edgeImg,1), size(edgeImg,2), filterCount);
    filters = cell(filterCount,1);
    for filtItr = 0:(filterCount-1)
        
       % If the filter does not exist, create it.
      if ~exist([currentFolder '/filters/filt' num2str(filtItr+1) '.png'], 'file') 
           theta = (-180/filterCount) * filtItr;
           gaborFilt = imrotate(refGaborFilt, theta, 'bilinear', 'crop');

           % Normalize the filter so there are no discrepancies between
           % them
           curConst = sum(sum(gaborFilt));
 %          gaborFilt = gaborFilt * (normConst/curConst);
           gaborFilt = gaborFilt - curConst/(numel(refGaborFilt));
           filters{filtItr+1} = gaborFilt;
           
           % Save filters
     %      imwrite(gaborFilt, [currentFolder '/filters/filt' num2str(filtItr+1) '.png']);
     %      imwrite(gaborFilt>gaborAreaMinResponse, [currentFolder '/filters/filt' num2str(filtItr+1) 'Mask.png']);
      else
   %       gaborFilt = double(imread([currentFolder '/filters/filt' num2str(filtItr+1) '.png']));
            load([currentFolder '/filters/filt' num2str(filtItr+1) '.mat']);
            filters{filtItr+1} = gaborFilt;
      end
       responseImg = conv2(double(edgeImg), gaborFilt, 'same');
       responseImg(responseImg<gaborFilterThr) = 0;
       testResp = @(x) testResponse(x);
       filteredImg = responseImg;
 %      filteredImg = colfilt(responseImg, [gaborFilterSize, gaborFilterSize], 'sliding', testResp);
       % Save response for future processing
       responseImgs(:,:,(filtItr+1)) = filteredImg;
    end
    
   %% Write smooth object boundaries based to an image on responseImgs.
   smoothedImg = getSmoothedImage(responseImgs, filters);
   
   %% Inhibit weak responses in vicinity of powerful peaks.
   peaks = find(responseImgs);
   [xInd, yInd, ~] = ind2sub(size(responseImgs), peaks);
   halfSize = floor(gaborFilterSize/2);
   
   for peakItr = 1:size(peaks,1)
      % If this peak has not yet been eliminated, go check nearby peaks
      if responseImgs(peaks(peakItr)) > 0 
          nearbyPeakIdx = xInd >= (xInd(peakItr) - halfSize) & xInd <= (xInd(peakItr) + halfSize) & ...
              yInd >= (yInd(peakItr) - halfSize) & yInd <= (yInd(peakItr) + halfSize);
          
          maxPeak = max(responseImgs(peaks(nearbyPeakIdx)));
          if responseImgs(peaks(peakItr)) > (maxPeak-0.00001)
             nearbyPeakIdx(peakItr) = 0;
             responseImgs(peaks(nearbyPeakIdx)) = 0; 
          else
             responseImgs(peaks(peakItr)) = 0;
          end
      end
   end
   
   %% Surpress weak responses.
   gaborThr = gaborAreaMinResponse * max(unique(responseImgs));
   responseImgs(responseImgs<gaborThr) = 0;
   
   % Write the responses in the final image.
   responseImgs = double(responseImgs>0);
   for filtItr = 1:filterCount
      responseImgs(:,:,filtItr) = responseImgs(:,:,filtItr) .* filtItr;
   end
   responseImg = sum(responseImgs,3);
   
   %% Out of this response image, we will create the nodes and output them.
   finalNodeIdx = find(responseImg);
   nodes = cell(numel(finalNodeIdx), 2);
   for nodeItr = 1:numel(finalNodeIdx)
       [centerX, centerY] = ind2sub(size(responseImg), finalNodeIdx(nodeItr));
       nodes(nodeItr,:) = {responseImg(finalNodeIdx(nodeItr)), [centerX, centerY]};
   end
end


%> Name: testResponse
%>
%> Description: Given a neighborhood, check if seed pixel has maximal
%> value. If it does, its value is returned. If not, 0 returned.
%>
%> @param nhood Column-wise neighborhood information (for optimization)
%>
%> @retval values Column-wise tested response.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 03.12.2013
function [values] = testResponse(nhood)
    [maxVal, maxIdx] = max(nhood, [], 1);
    areEqual = maxVal >= ((nhood(ceil(size(nhood,1)/2),:))-0.00001) & maxVal > 0;
    values = (maxVal .* areEqual)';
end



function [smoothedImg] = getSmoothedImage(responseImgs, filters)
    smoothedImg = zeros(size(responseImgs,1), size(responseImgs,2));
    halfSize = floor(size(filters{1},1)/2);
    responseImgs([1:(halfSize+1), (end-halfSize):end], :, :) = 0;
    responseImgs(:, [1:(halfSize+1), (end-halfSize):end], :) = 0;
    responseImgs = responseImgs/max(max(max(responseImgs)));
    for filterItr = 1:size(filters,1)
        [xInd, yInd, val] = find(responseImgs(:,:,filterItr));
        for peakItr = 1:size(peaks,1)
            smoothedImg((xInd(peakItr)-halfSize):(xInd(peakItr)+halfSize), ...
                (yInd(peakItr)-halfSize):(yInd(peakItr)+halfSize)) = val(peakItr) * filters{filterItr};
        end
 %       smoothedImg = smoothedImg + conv2(responseImgs(:,:,filterItr), filters{filterItr});
    end
    smoothedImg = smoothedImg/max(max(smoothedImg));
end

