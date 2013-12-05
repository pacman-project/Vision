% this is to test a new line discretization function

root = 'D:\3D\LibHoP3D\';
addpath([root,'utility']);
addpath([root,'settings']);

dxKernel4 = [1,-8,0,8,-1]; % consistancy order 4 
dxKernel4 = dxKernel4 / 12;
kernelSize = 3;
sigma = 0.6;
nClusters = 9;
thresh = 95;

% ------------ parameters for line discretization---------------------
% min [error + wOverlap * overlap - wCoverage*coverage ]
wCoverage = 1;
wOverlap = 2;
% ------------ parameters for line discretization---------------------

kernelW = load('settings/kernelW.mat');
kernelW = kernelW.kernelW;    % this is the kernel for lineDiscretization and lineDiscretizationY

% cluster sizes in percent (if we have nClusters = 9)
if nClusters == 9
    clusterSizes = load('settings/cluster1Sizes.mat'); 
    clusterSizes = clusterSizes.clusterSizes;
    cluster1SizesPercent = clusterSizes* 0.01;
end

combs = load('settings/combs12.mat');
combs = combs.combs;

[cluster1Centers, cluster1Length] = defineCluster1Centers(nClusters, cluster1SizesPercent, thresh);

I = imread('D:\3D\Input Data\Images for categorization\1view_1Scale\images\1_25_1_1_1_1.png');
I = double(I);

[r,c] = size(I);
mask = zeros(r,c);
mask(I > 0) = 1;

% check if the mask is empty or not
maxM = max(max(mask));
if maxM == 0 
    return;
end   

[I, mask] = trimImage(I, mask); % trim the image

% compute filter responces
hx = dxKernel4;
hy = hx'; 
smallestLine = 10; % the smallest line used for learning

h = fspecial('gaussian', kernelSize, sigma);
IG = imfilter(I, h); % gaussian                 % done

% now computing Gaussian derivatives
Ix = imfilter(IG, hx);
Ix = Ix.*mask; % to avoid high on the boundary

% perform line discretization (line after line)
[r, c] = size(Ix);
marks = zeros(r,c); % here we collect marks


for k = 2:2:100 % for every line we perform discretization
    ind = find(mask(k,:));
    if length(ind) < smallestLine
        marks(k,:) = zeros(1, c); % line of zeros in this case
    else
        fx = Ix(k,ind(1):ind(end));
        output = lineDiscretizationOptLayer1(fx, wCoverage, wOverlap, nClusters, cluster1Centers, cluster1Length, thresh, combs);
        marks(k,ind(1):ind(end)) = output;
    end
    k
end
imtool(marks, [0, nClusters]);  %   GIVES NICE PICTURES 
marks = marks .* mask;




































