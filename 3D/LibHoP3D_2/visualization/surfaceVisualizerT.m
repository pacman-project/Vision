% this function is to visualize any 3 layer elements.

% positions: [x1,y1,z1; x2,y2,z2; ...] may be real numbers far from margin
% elements (35,44,51,...) - index of the corresponding second layer part
% nCluster - integer
% thresh - real number (refers to thresh form the 1st layer)
% fieldSize is a vector [sizeX, SizeY, sizeZ]
% depthStep - real number (presumably 95/5)

function [out] = surfaceVisualizerT(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep, faceColor, isFlip, transparency)

    if nargin < 9
        transparency = 1.0;
        edgeTransparency = 1.0;
    else
        edgeTransparency = 0.0;
    end
        
    % set default parameters here
    if (nargin == 6)
        faceColor = [0 0 1];
        isFlip = false;
    elseif (nargin == 7)
        isFlip = false;
    end
    
    n2Clusters = nClusters^2;
    emptyID = n2Clusters + 1;
        

    len = length(elements);
    lenV = len * 4; % two triangles for each part
    lenF = len * 2;
    
    V = zeros(lenV, 3);
    F = zeros(lenF, 3);
    
    offset = 1;
    
    if fieldSize(1) > fieldSize(2)
        maxDim = fieldSize(1) + 2*offset;
    else
        maxDim = fieldSize(2) + 2*offset;
    end
    
    step = 2.5;
    
    for i = 1:len
        if elements(i) == emptyID
            
            V((i-1)*4 + 1, :) = [1,1,1];  % vertex
            V((i-1)*4 + 2, :) = [1,1,1];  % vertex  
            V((i-1)*4 + 3, :) = [1,1,1];  % vertex  
            V((i-1)*4 + 4, :) = [1,1,1];  % vertex 
            
        else
            
            [clusterX, clusterY] = compute2derivatives(elements(i), nClusters);  
            dx = cluster1Centres(clusterX);
            dy = cluster1Centres(clusterY);

            curCentre = positions(i, :);
            curCentre(3) = curCentre(3) * depthStep;

            V((i-1)*4 + 1, :) = [curCentre(1) - step, curCentre(2) - step, curCentre(3) - step*dx - step*dy];  % vertex
            V((i-1)*4 + 2, :) = [curCentre(1) - step, curCentre(2) + step, curCentre(3) - step*dx + step*dy];  % vertex  
            V((i-1)*4 + 3, :) = [curCentre(1) + step, curCentre(2) + step, curCentre(3) + step*dx + step*dy];  % vertex  
            V((i-1)*4 + 4, :) = [curCentre(1) + step, curCentre(2) - step, curCentre(3) + step*dx - step*dy];  % vertex  
        end
        
        F((i-1)*2 + 1, :) = [(i-1)*4 + 1, (i-1)*4 + 2, (i-1)*4 + 3]; % faces
        F((i-1)*2 + 2, :) = [(i-1)*4 + 1, (i-1)*4 + 3, (i-1)*4 + 4];

    end
    
    % make some offsets  but why do we need them
%     
%     V(:, 1:2) = V(:, 1:2) + 1;

    trisurf(F, V(:,1),V(:,2),V(:,3),'FaceColor',faceColor, 'EdgeColor', faceColor, 'FaceAlpha', transparency, 'EdgeAlpha', edgeTransparency);
    xlabel('x')
    ylabel('y')
    axis([1, maxDim, 1, maxDim, 1, ceil(depthStep * fieldSize(3)), 50, 60]); % ceil(1 + scale * rangeZ), ceil(rangeZ - scale * rangeZ)
    if isFlip
        set(gca,'zdir','reverse')
    end
    
    
    out = true;
end
