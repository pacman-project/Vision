% this function performs actual reconstruction of the part to surface
% descriptors

% we assume the following levels of downsampling:
% 1) only central pixels of each part are taken

% 00000
% 00000
% 00100
% 00000
% 00000

% 2) only four pixels are taken

% 00000
% 01010
% 00000
% 01010
% 00000

% 3) many pixels are taken into account pixels are taken

% 10101
% 01010
% 10101
% 01010
% 10101

function [X, d] = partToSurfaceDescriptor(elements, positions, nClusters, cluster1Centres, downsamplingScheme)
    
    lenE = length(elements);

    if downsamplingScheme == 1 
        X = zeros(1, lenE);
    elseif downsamplingScheme == 2 
        X = zeros(1, 4 * lenE);
    elseif downsamplingScheme == 3 
        X = zeros(1, 13 * lenE);
    end
    
    for i = 1 : lenE
        
        curEl = elements(i);
        
        [clusterX, clusterY] = compute2derivatives(elements(i), nClusters);
        dx = cluster1Centres(clusterX);
        dy = cluster1Centres(clusterY);
        
        if downsamplingScheme == 1
            X(i) = positions(i, 3);
        elseif downsamplingScheme == 2
            
            X((i-1)*4 + 1) = positions(i, 3) - dx - dy;
            X((i-1)*4 + 2) = positions(i, 3) + dx - dy;
            X((i-1)*4 + 3) = positions(i, 3) + dx + dy;
            X((i-1)*4 + 4) = positions(i, 3) - dx + dy;
            
        elseif downsamplingScheme == 3
            
            X((i-1)*13 + 1) = positions(i, 3) - 2*dx - 2*dy;
            X((i-1)*13 + 2) = positions(i, 3) + 0    - 2*dy;
            X((i-1)*13 + 3) = positions(i, 3) + 2*dx - 2*dy;
            X((i-1)*13 + 4) = positions(i, 3) - dx - dy;
            X((i-1)*13 + 5) = positions(i, 3) + dx - dy;
            X((i-1)*13 + 6) = positions(i, 3) - 2*dx ;
            X((i-1)*13 + 7) = positions(i, 3) - 0 ;
            X((i-1)*13 + 8) = positions(i, 3) - 2*dx ;
            X((i-1)*13 + 9) = positions(i, 3) - dx + dy;
            X((i-1)*13 + 10)= positions(i, 3) + dx + dy;
            X((i-1)*13 + 11)= positions(i, 3) - 2*dx + 2*dy;
            X((i-1)*13 + 12)= positions(i, 3) - 0 + 2*dy;
            X((i-1)*13 + 13)= positions(i, 3) + 2*dx + 2*dy;
        end
    
    end
    
    d = length(X);

end

























