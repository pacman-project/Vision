% this is to test a new line discretization function

dxKernel4 = [1,-8,0,8,-1]; % consistancy order 4 
dxKernel4 = dxKernel4 / 12;
kernelSize = 3;
sigma = 0.6;

I = imread('D:\3D\Input Data\Images for categorization\1view_1Scale\images\1_38_1_1_1_1.png');
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


for k = 2:2:r % for every second line we perform discretization
    ind = find(mask(k,:));
    if length(ind) < smallestLine
        marks(k,:) = zeros(1, c); % line of zeros in this case
    else
        fx = Ix(k,ind(1):ind(end));
        strlen = length(fx); % strlen must be greater than smallestLine
        output = lineDiscretization(fx, strlen, nClusters, clusterCenters, clusterSize, thresh);
        outline = zeros(1,c);
        marks(k,ind(1):ind(end)) = output;
    end
end
imtool(marks, [1,nClusters]);  %   GIVES NICE PICTURES 
marks = marks .* mask;



