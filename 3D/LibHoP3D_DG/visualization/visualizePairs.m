% this function is for visualization of pairs
% positions should be in the following format: [x,y,z, Nx, Ny, Nz]

function is_ok = visualizePairs(positions)

    if nargin == 1
        transparency = 1.0;
        edgeTransparency = 1.0;
        faceColor = [0 0 1];
        isFlip = false;
    end
    len = size(position, 1);
        
    for i = 1:len
            
            dx = positions(i, );%atan(cluster1Centres(clusterX)*pi/180);
            dy = positions();%atan(cluster1Centres(clusterY)*pi/180);

            curCentre = positions(i, :);
%             curCentre(3) = curCentre(3) * depthStep;

            V((i-1)*4 + 1, :) = [curCentre(1) - step, curCentre(2) - step, curCentre(3) - step*dx - step*dy];  % vertex
            V((i-1)*4 + 2, :) = [curCentre(1) - step, curCentre(2) + step, curCentre(3) - step*dx + step*dy];  % vertex  
            V((i-1)*4 + 3, :) = [curCentre(1) + step, curCentre(2) + step, curCentre(3) + step*dx + step*dy];  % vertex  
            V((i-1)*4 + 4, :) = [curCentre(1) + step, curCentre(2) - step, curCentre(3) + step*dx - step*dy];  % vertex  
        end
        
        F((i-1)*2 + 1, :) = [(i-1)*4 + 1, (i-1)*4 + 2, (i-1)*4 + 3]; % faces
        F((i-1)*2 + 2, :) = [(i-1)*4 + 1, (i-1)*4 + 3, (i-1)*4 + 4];

    end
    
    trisurf(F, V(:,1),V(:,2),V(:,3),'FaceColor',faceColor, 'EdgeColor', faceColor, 'FaceAlpha', transparency, 'EdgeAlpha', edgeTransparency);
    xlabel('x')
    ylabel('y') 
    axis([0, maxDim, 0, maxDim, 0, maxDim, 50, 60]); % ceil(1 + scale * rangeZ), ceil(rangeZ - scale * rangeZ)
    if isFlip
        set(gca,'zdir','reverse')
    end
    

end

