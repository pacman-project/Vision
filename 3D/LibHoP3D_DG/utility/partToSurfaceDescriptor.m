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

% in the emptyIndicator empty cells are marked with ONES

function [X, emptyIndicator] = partToSurfaceDescriptor(elements, positions, nClusters, cluster1Centres, downsamplingScheme)
    
    lenE = length(elements);

    if downsamplingScheme == 1 
        X = zeros(1, lenE);
        emptyIndicator = zeros(1, lenE);
    elseif downsamplingScheme == 2 
        X = zeros(1, 4 * lenE);
        emptyIndicator = zeros(1, 4 * lenE);
    elseif downsamplingScheme == 3 
        X = zeros(1, 13 * lenE);
        emptyIndicator = zeros(1, 13 * lenE);
    end
    
    for i = 1 : lenE
        
        curEl = elements(i);
        
        if curEl == nClusters^2 + 1 % this is a detection of an empty cell
            if downsamplingScheme == 1
                X(i) = 0;
                emptyIndicator(i) = 1;
            elseif downsamplingScheme == 2
                X((i-1)*4 + 1:(i-1)*4 + 4) = zeros(1, 4);
                emptyIndicator((i-1)*4 + 1:(i-1)*4 + 4) = ones(1, 4);
            elseif downsamplingScheme == 3
                X((i-1)*13 + 1:(i-1)*13 + 13) = zeros(1, 13);
                emptyIndicator((i-1)*13 + 1:(i-1)*13 + 13) = ones(1, 13);
            end
        else
            [clusterX, clusterY] = compute2derivatives(elements(i), nClusters);
            
            % dx and dy should be slopes  (dz/dx) and (dz/dy)
            dx = atan(cluster1Centres(clusterX)*pi/180);
            dy = atan(cluster1Centres(clusterY)*pi/180);

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
    
    end

end

























