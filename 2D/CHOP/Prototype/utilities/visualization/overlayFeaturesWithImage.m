function [img] = overlayFeaturesWithImage(level1Nodes, img, filters)
    filterSize = size(filters{1});
    filterWeight = 0.8;
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
    filters = cellfun(@(x) uint8(255 * ((x - minVal)/(maxVal - minVal))), filters, 'UniformOutput', false);
    
    for nodeItr = 1:size(level1Nodes,1)
        lowX = level1Nodes(nodeItr,2) - lowHalf(1);
        highX = level1Nodes(nodeItr,2) + highHalf(1);
        lowY = level1Nodes(nodeItr,3) - lowHalf(2);
        highY = level1Nodes(nodeItr,3) + highHalf(2);
        
        if lowX < 1 || lowY < 1 || highX>size(img,1) || highY>size(img,2)
            continue;
        end
        tempImg = img;
        tempImg(lowX:highX, lowY:highY, :) = uint8((imageWeight * double(orgImg(lowX:highX, lowY:highY, :)) + ...
             filterWeight * double(filters{level1Nodes(nodeItr,1)})));
        img = max(tempImg, img);
    end
end

