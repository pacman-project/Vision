% in this script we analise the distribution in the first layer of
% hierarchy and select the base element

% interpretation of depths: range from 0 to 20 in terms of global coords
% [allDXs, allDYs, lenX, LenY] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel, histSize);

function [allDXs, lenX] = computeFirstLayerDistributionMesh(list_depth)

    disp('Learning statistics of the first layer. Iterations started...');
    lenF = length(list_depth);
    allDXs = [];

    parfor k = 1 : lenF
        [V, F, N] = meshRead(filename);
        

        [~, Ix, ~, mask, ~, ~, is_successfull] = preliminaryProcessingMesh(V, F, N, options);

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







