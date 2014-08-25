I = imread('F:\Washington Cropped\rgbd-dataset\apple\apple_1\apple_1_1_1_depthcrop.png');
mask = imread('F:\Washington Cropped\rgbd-dataset\apple\apple_1\apple_1_1_1_maskcrop.png');
I = imresize(I, scale);
mask = imresize(mask, scale);
scale = 2.5;

% gaussian
kernel =fspecial('gaussian', 25, 6.0);
I = double(I);
mask = double(mask);
Ig = imfilter(I, kernel, 'replicate', 'same');


% another mask
[r,c] = size(I);
mask1 = ones(r,c);
mask1(I < 600 & mask == 1) = 0;
mask1(I > 800 & mask == 1) = 0;

se = strel('disk',2);

imtool(mask);
mask(mask1 == 0) = 0;
mask = imerode(mask, se);
imtool(mask);

% filtering of the image under guidance of mask

%   GUIDEDFILTER   O(1) time implementation of guided filter.
%
%   - guidance image: I (should be a gray-scale/single channel image)
%   - filtering input image: p (should be a gray-scale/single channel image)
%   - local window radius: r
%   - regularization parameter: eps

r = 8;
eps = 0.05^2;

Iguided = guidedfilter(mask, I, r, eps);

% measure errors:
imtool(Ig, [min(min(Ig)), max(max(Ig))]);
imtool(Iguided, [min(min(Iguided)), max(max(Iguided))]);

Ig(mask == 0) = 0;
Iguided(mask == 0) = 0;

errG = sum(sum(abs(I - Ig)))/r/c 
errGu = sum(sum(abs(I - Iguided)))/r/c




