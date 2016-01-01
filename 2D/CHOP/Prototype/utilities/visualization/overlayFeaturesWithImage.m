%> Name: overlayFeaturesWithImage
%>
%> Description:Given a set of level 1 nodes, the actual image, and filter
%> images, this function overlays the level 1 filters with the original
%> image. The image's weight is 1-filterWeight, while the filters are put 
%> into the original image with weight filterWeight. 
%>
%> @param level1Nodes Level 1 nodes, an Nx5 array that has the following
%> format: [label, posX, posY, levelItr, imageId; ...].
%> @param img Original image.
%> @filters The cell vector containing filter images.
%> filterWeight Between 0-1. It specifies the filter weights when overlaid
%> with the original image, with the original image contributing to the final
%> image with a factor of 1-filterWeight.
%>
%> @retval img The overlay image that visualizes both parts.
%> 
%> Author: Rusen
%>
%> Updates
%> Ver 1.0 on 19.11.2015
function [img] = overlayFeaturesWithImage(level1Nodes, img, filters, filterWeight)
    filterSize = size(filters{1});
    imageWeight = 1 - filterWeight;
    centerPoint = round(filterSize/2);
    lowHalf = centerPoint-1;
    highHalf = filterSize - centerPoint;
    orgImg = img;
    if max(max(max(img))) == 1
       img = uint8(round(double(img) * 255));
    end
    img = uint8(round(imageWeight * double(img)));
 
    minVal = min(cellfun(@(x) min(min(min(x))), filters));
    maxVal = max(cellfun(@(x) max(max(max(x))), filters));
    if minVal < 0 || maxVal <= 1
          filters = cellfun(@(x) uint8(255 * ((double(x) - minVal)/(maxVal - minVal))), filters, 'UniformOutput', false);
    end
    imgCounts = zeros(size(img,1), size(img,2));
    tempImg = zeros(size(img));
    for nodeItr = 1:size(level1Nodes,1)
        lowX = level1Nodes(nodeItr,2) - lowHalf(1);
        highX = level1Nodes(nodeItr,2) + highHalf(1);
        lowY = level1Nodes(nodeItr,3) - lowHalf(2);
        highY = level1Nodes(nodeItr,3) + highHalf(2);
        
        if lowX < 1 || lowY < 1 || highX>size(img,1) || highY>size(img,2)
            continue;
        end
        tempImg(lowX:highX, lowY:highY, :) = tempImg(lowX:highX, lowY:highY, :) + double(filters{level1Nodes(nodeItr,1)});
        imgCounts(lowX:highX, lowY:highY) = imgCounts(lowX:highX, lowY:highY) + ones(filterSize(1), filterSize(2));
    end
    imgCounts(imgCounts == 0) = 1;
    for bandItr = 1:size(tempImg,3)
         tempBandImg = tempImg(:,:,bandItr);
         tempBandImg = tempBandImg ./ imgCounts;
         tempImg(:,:,bandItr) = tempBandImg;
    end
    
    img = uint8(round(tempImg * filterWeight + double(orgImg) * imageWeight));
end

