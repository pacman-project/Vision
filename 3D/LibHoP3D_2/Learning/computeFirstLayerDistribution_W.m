% in this script we analise the distribution in the first layer of
% hierarchy and select the base element

% interpretation of depths: range from 0 to 20 in terms of global coords
% [allDXs, allDYs, lenX, LenY] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel, histSize);

function [allDXs, lenX] = computeFirstLayerDistribution_W(list_depth, list_mask, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, ...
                               is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2)

    disp('Learning statistics of the first layer. Iterations started...');
    lenF = length(list_depth);
    allDXs = [];
    
    isTrim = true;
    isY = false;
    isX = true;

    parfor k = 1 : lenF
        I = imread(list_depth{k});
        mask = imread(list_mask{k});
        
        I = I(:,:,1);             
        [~, Ix, ~, mask, ~, ~, is_successfull] = preliminaryProcessing(I, mask, isErrosion, discSize, isX, isY, isTrim, dxKernel,...
                        sigmaKernelSize, sigma, is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2);

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







