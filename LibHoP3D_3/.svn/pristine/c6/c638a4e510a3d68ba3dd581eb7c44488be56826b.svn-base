% this function is to visualize any 3 layer elements.

% positions: [x1,y1,z1; x2,y2,z2; ...] may be real numbers far from margin
% elements (35,44,51,...)
% nCluster - integer
% thresh - real number (refers to thresh form the 1st layer)
% fieldSize is a vector [sizeX, SizeY, sizeZ]
% depthStep - real number (presumably 95/5)

function [out] = surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep)

    len = length(elements);
    
    if fieldSize(1) > fieldSize(2)
        maxDim = fieldSize(1)+5;
    else
        maxDim = fieldSize(2)+5;
    end
    
    Xbig = NaN(maxDim, maxDim);
    Ybig = NaN(maxDim, maxDim);
    Zbig = NaN(maxDim, maxDim);
 
%     Xbig = NaN(fieldSize(1), fieldSize(2));
%     Ybig = NaN(fieldSize(1), fieldSize(2));
%     Zbig = NaN(fieldSize(1), fieldSize(2));
     
    [Xmesh, Ymesh] = meshgrid(1:1:maxDim);
    
    for i = 1:len    
        [clusterX, clusterY] = compute2derivatives(elements(i), nClusters);  

        [X, Y] = meshgrid(1:1:5);
        dx = cluster1Centres(clusterX);
        dy = cluster1Centres(clusterY);
        
        Z = dx * (X) + dy * Y;
        Z = Z - Z(3,3);
        Z = Z + positions(i,3)*depthStep;
        
        % then put it to the right position in the file
        curShiftX = round(positions(i,1));
        curShiftY = round(positions(i,2));
        
        patchSize = 2;
        if(curShiftX>patchSize)&&(curShiftY>patchSize)&&(curShiftX<maxDim-patchSize)&&(curShiftY<maxDim-patchSize)
        
            Xbig(curShiftY-2:curShiftY+2, curShiftX-2:curShiftX+2) = Xmesh(curShiftY-2:curShiftY+2, curShiftX-2:curShiftX+2);
            Ybig(curShiftY-2:curShiftY+2, curShiftX-2:curShiftX+2) = Ymesh(curShiftY-2:curShiftY+2, curShiftX-2:curShiftX+2);
            Zbig(curShiftY-2:curShiftY+2, curShiftX-2:curShiftX+2) = Z;
            
        end;
    end

    % mesh(Xbig,Ybig,Zbig, 'EdgeColor','blue');
    surf(Xbig,Ybig,Zbig);
    xlabel('x')
    ylabel('y')
    axis([1, fieldSize(1), 1, fieldSize(2), 1, ceil(depthStep * fieldSize(3)), 50, 60]); % ok. done
    
    out = true;
end
