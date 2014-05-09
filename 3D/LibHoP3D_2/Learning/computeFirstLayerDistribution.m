% in this script we analise the distribution in the first layer of
% hierarchy and select the base element

% interpretation of depths: range from 0 to 20 in terms of global coords
% [allDXs, allDYs, lenX, LenY] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel, histSize);

function [allDXs, lenX] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, is_downsampling, dowsample_rate)


    disp('Learning statistics of the first layer. Iterations started...');
    len = length(list_depth);
    allDXs = [];
    
    isTrim = true;
    isY = false;
    isX = true;

    parfor k = 1 : len
        I = imread(list_depth{k});
        
        if is_downsampling
            I = imresize(I, dowsample_rate);
        end
      
        [~, Ix, ~, mask, ~, ~, is_successfull] = preliminaryProcessing(I, [], isErrosion, discSize, isX, isY, isTrim, dxKernel, sigmaKernelSize, sigma)

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







