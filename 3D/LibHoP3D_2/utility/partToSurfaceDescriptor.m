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

% 3) all four pixels are taken

% 11111
% 11111
% 11111
% 11111
% 11111

function [X, d] = partToSurfaceDescriptor(elements, positions, nClusters, cluster1Centres, downsamplingScheme)
    
    lenE = length(elements);

    if downsamplingScheme == 1 
        X = zeros(1, lenE);
    elseif downsamplingScheme == 2 
        X = zeros(1, 4 * lenE);
    elseif downsamplingScheme == 3 
        X = zeros(1, 15 * lenE);
    end
    
    for i = 1 : length(elements)
        
        curEl = elements(i);
        
        if downsamplingScheme == 1
            X(i) = positions(i, 3);
        elseif downsamplingScheme == 2
            [clusterX, clusterY] = compute2derivatives(elements(i), nClusters);
            dx = cluster1Centres(clusterX);
            dy = cluster1Centres(clusterY);
            
            X((i-1)*4 + 1) = positions(i, 3) + dx + dy;
            X((i-1)*4 + 2) = positions(i, 3) + dx - dy;
            X((i-1)*4 + 3) = positions(i, 3) - dx - dy;
            X((i-1)*4 + 4) = positions(i, 3) - dx + dy;
            
        elseif downsamplingScheme == 3
        
            disp('ERROR: this downsampling schema is not defined!!!')
        end
    
    end
    
    d = length(X);

end

























