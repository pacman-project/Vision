dataSetNumber = 2; 

I = imread('D:\Input Data\Washington\Wash-rgbd-dataset\ball\ball_2\ball_2_1_1_depthcrop.png');

I = I(:,:,1);
I = double(I);
[r,c] = size(I);
imtool(I, [min(min(I)), max(max(I))]);
 
 % take the gaussian deriveatives
 
[dxKernel, combs, largestLine, sigma, sigmaKernelSize, isErrosion, discRadius, is_guided, r_guided, eps, ...
is_mask_extended, maxExtThresh1, maxExtThresh2] = loadFilteringParameters(dataSetNumber);
sigma = 2.0;
sigmaKernelSize = 9;

h = fspecial('gaussian', sigmaKernelSize, sigma);
IG = imfilter(I, h, 'replicate', 'same'); % gaussian

hx = dxKernel; 
Ix = imfilter(IG, hx, 'replicate', 'same');
 
I(I < 690) = NaN;
I(I > 900) = NaN;
I(Ix < -2) = NaN;
I(Ix > 2)  = NaN;

% extract the table

leftRightKernel = [];



 
 imtool(I, [min(min(I)), max(max(I))]);
