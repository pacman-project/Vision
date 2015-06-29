% this script is to test the structure tensor on depth images

is2D = false;
dataSetNumber = 3;

[filtOptions, zeroThresh] = loadFilteringParameters(dataSetNumber);
options = GetOptions(filtOptions);

I = imread('D:\Input Data\VladislavSTD\Vladislav_STD\depthEdgesCorners\99_1_1_1_2_3.png');

load('settings/calibration3.mat');
I = double(I(:,:,1));
I= I * zScale;
[r,c, ~] = size(I);
mask = zeros(r,c);
mask(I > 0) = 1;

dx = [0.5, 0, -0.5];
dy = dx';
h = fspecial('gaussian', 5, 1);
I = imfilter(I, h);

Ix = imfilter(I, dx);
Iy = imfilter(I, dy);

Ix(mask == 0) = 0;
Iy(mask == 0) = 0;
I(mask == 0) = 0;

borderOffset = 10;

Normals = zeros(r,c,3);
% compute normals at each point
parfor x = borderOffset:c-borderOffset
    for y = borderOffset:r-borderOffset
        
        if mask(y,x) == 0
            continue;
        end
        
        % compute tangent vectors
        Tx = [1, 0, Ix(y,x)];
        Ty = [0, 1, Iy(y,x)];
        Tx = Tx/norm(Tx);
        Ty = Ty/norm(Ty);
        
        N = cross(Tx, Ty);
        for j = 1:3
            Normals(y,x,j) = N(j);
        end
        
    end
end

Nx = Normals(:,:,1);
Ny = Normals(:,:,2);
Nz = Normals(:,:,3);
imtool(Normals);

Iout = zeros(r,c,3);

parfor x = borderOffset:c-borderOffset
    for y = borderOffset:r-borderOffset
        
        if mask(y,x) == 0
            continue;
        end
        
        [indsXOut, indsYOut, ~] = computeNeighbors(I, x, y, borderOffset-1, is2D, filtOptions.minDepth, filtOptions.maxDepth);
        if length(indsXOut) < 45
            continue;
        end
        
        inds = sub2ind([r,c], indsYOut, indsXOut);
        Ns_x = Nx(inds)';
        Ns_y = Ny(inds)';
        Ns_z = Nz(inds)';
        Ns = [Ns_x, Ns_y, Ns_z];
        
        [T, V, D] = structureTensor3(Ns);
        b = diag(D);
        
        for j =1:3
            Iout(y,x,j) = b(j);
        end
 
    end
end

% here we estimate the local frame of reference for the edge-corner parts



imtool(Iout);




