% this function is to visualize any 3 layer elements.

% positions: [x1,y1,z1; x2,y2,z2; ...] may be real numbers far from margin
% elements (35,44,51,...)
% nCluster - integer
% thresh - real number (refers to thresh form the 1st layer)
% fieldSize is a vector [sizeX, SizeY, sizeZ]
% depthStep - real number (presumably 95/5)

function [out] = surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep, faceColor)

    if (nargin == 6)
        faceColor = [0 0 1];
    end

    len = length(elements);
    
    offset = 1;
    
    if fieldSize(1) > fieldSize(2)
        maxDim = fieldSize(1) + 2*offset;
    else
        maxDim = fieldSize(2) + 2*offset;
    end
    
    [Xmesh, Ymesh] = meshgrid(1:0.5:maxDim);
    
    meshSize = size(Xmesh, 1);
    Xbig = NaN(meshSize, meshSize);
    Ybig = NaN(meshSize, meshSize);
    Zbig = NaN(meshSize, meshSize);
    
    elementSize = 6;
    halfElSize = elementSize - 1;
    elCentre = elementSize;
    
    for i = 1:len    
        [clusterX, clusterY] = compute2derivatives(elements(i), nClusters);  

        [X, Y] = meshgrid(1:0.5:elementSize);
        dx = cluster1Centres(clusterX);
        dy = cluster1Centres(clusterY);
        
        Z = dx * (X) + dy * Y;
        Z = Z - Z(elCentre,elCentre);
        Z = Z + positions(i,3)*depthStep;
        
        % then put it to the right position in the file
        curShiftX = round(positions(i,1)) * 2 + offset;
        curShiftY = round(positions(i,2)) * 2 + offset;
        
        Xbig(curShiftY-halfElSize:curShiftY+halfElSize, curShiftX-halfElSize:curShiftX+halfElSize) = Xmesh(curShiftY-halfElSize:curShiftY+halfElSize, curShiftX-halfElSize:curShiftX+halfElSize);
        Ybig(curShiftY-halfElSize:curShiftY+halfElSize, curShiftX-halfElSize:curShiftX+halfElSize) = Ymesh(curShiftY-halfElSize:curShiftY+halfElSize, curShiftX-halfElSize:curShiftX+halfElSize);
        Zbig(curShiftY-halfElSize:curShiftY+halfElSize, curShiftX-halfElSize:curShiftX+halfElSize) = Z;

    end

    surf(Xbig,Ybig,Zbig, 'EdgeColor', 'none', 'FaceColor',faceColor,'FaceAlpha',0.8);
    xlabel('x')
    ylabel('y')
    axis([1, maxDim, 1, maxDim, 1, ceil(depthStep * fieldSize(3)), 50, 60]); % ceil(1 + scale * rangeZ), ceil(rangeZ - scale * rangeZ)
    
    out = true;
end
