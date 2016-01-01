% this is a function that converts a line of format
% lineIn:  central, leftRO, relDepthL, rightRO, relDepthR, partID_left, partID_right

 % lineout: to the format central, left, relDepthL, right, relDepthR

function [lineout, ConversionTable] = DiscretizeOrientations(lineIn, layerID, curTS, nClusters)

    if layerID == 3
        centralPart = compute2elementIndex(ceil(nClusters/2), ceil(nClusters/2), nClusters);
        lineout = zeros(curTS, 5);
        lineout(:, 1) = ones(curTS, 1) * centralPart;
        lineout(:, 2:5) = lineIn(:, 2:5);
        ConversionTable = [];
    end
    
end


