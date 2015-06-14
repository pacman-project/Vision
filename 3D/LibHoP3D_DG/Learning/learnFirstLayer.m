% this function is purposed to learn the parameters of the first layer
% applicable for the Washington dataset only

function [cluster1Centres, cluster1Bounds, thresh] = learnFirstLayer(list_depth, list_mask, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, ...
                      nClusters, quantilesFirst, is_guided, r_guided, eps, ...
                      is_mask_extended, maxExtThresh1, maxExtThresh2)
    
    if ~((nClusters == 9)||(nClusters == 7))
        disp('Error. Please define quantiles for this number of clusters');
    end
    
    allDXs = zeros(10^7,1); % here we collect all the values of dx
    
    [allDXs, lenX] = computeFirstLayerDistribution(list_depth, list_mask, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize,...
                          is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2);
 
    allDXs = allDXs(1:lenX);
    
    if nClusters == 9 || nClusters == 7
        [bounds] = quantile(allDXs, quantilesFirst);
    end
    
    %% Here we process statistics defining cluster centres and lengths
    
    cluster1Bounds = [bounds, fliplr(abs(bounds))]; % to make it symmetric (if for whatever reason it is not)
    cluster1Centres = zeros(1, nClusters);
    for i = 1:nClusters
        cluster1Centres(i) = (cluster1Bounds(i) + cluster1Bounds(i+1)) / 2;
    end
    
%   cluster1Centres(1) = cluster1Bounds(2) - 0.01;
%   cluster1Centres(nClusters) = cluster1Bounds(nClusters) + 0.01;
    thresh = -cluster1Centres(1);
    
end
