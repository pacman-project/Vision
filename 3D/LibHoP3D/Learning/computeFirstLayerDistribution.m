% in this script we analise the distribution in the first layer of
% hierarchy and select the base element

% interpretation of depths: range from 0 to 20 in terms of global coords
% [allDXs, allDYs, lenX, LenY] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel, histSize);

function [allDXs, lenX] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel)

    disp('Learning statistics of the first layer. Iterations started...');
    len = length(list_depth);
    allDXs = [];
    
    isTrim = true;
    isErrosion = true;
    discSize = 3;
    isY = false;

    parfor k = 1 : len
        I = imread(list_depth{k});
      
        [I, Ix, Iy, mask, r, c, is_successfull] = preliminaryProcessing(I, isErrosion, discSize, isY, isTrim, dxKernel, sigmaKernelSize, sigma)

%         [r,c] = size(I);
%         I = double(I);
%         mask = zeros(r,c);
%         mask(I > 0) = 1;
% 
%         hx = dxKernel;
%         % hy = hx';       
%         h = fspecial('gaussian', sigmaKernelSize, sigma);
% 
%         % now computing Gaussian derivatives
%         IG = imfilter(I, h);
%         Ix = imfilter(IG, hx);
%         Ix = Ix.*mask; % to avoid high values on the boundary
% %       Iy = imfilter(IG, hy); % derivative in y direction
% %       Iy = Iy.*mask;
% 
%         % perform the errosion of the mask, with circle of size 2*sigma
%         se = strel('disk', 3);
%         mask = imerode(mask, se);
% 
%         Ix = Ix .* mask;

        if is_successfull
            Ia = Ix(mask == 1);
            allDXs = [allDXs; Ia];
        end

        if mod(k,20) == 0
            k % this is to trace the progress
        end
    end
    lenX = length(allDXs);
end






