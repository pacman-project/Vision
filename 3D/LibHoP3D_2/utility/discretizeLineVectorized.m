
% this function convert string of filter responces to the string of
% clucters and errors


function [nearestClusters, alternativeClusters] = discretizeLineVectorized(fx, nClusters, cluster1Centres, cluster1Bounds)
    
    % initialization
    alternativeClusters = zeros(size(fx));
    curCentres = zeros(size(fx));
    
    % compute nearestClusters and alternativeClusters
    nearestClusters = imquantize(fx, cluster1Bounds) - 1;
    nearestClusters(nearestClusters == nClusters + 1) = 0; % above the largest thersh
   
    % alternative clusters
    inds = nearestClusters > 0; % assigned to any bins
    curCentres(inds) = cluster1Centres(nearestClusters(inds));
    Shift = fx - curCentres;
    
    
    alternativeClusters(inds) = nearestClusters(inds) + sign(Shift(inds));
    alternativeClusters(alternativeClusters == nClusters + 1) = 0;
    
    
    alternativeClusters(fx < cluster1Centres(1)) = 1;
    alternativeClusters(fx > cluster1Centres(end)) = nClusters;
    
    alternativeClusters(fx < cluster1Bounds(1)) = 0;
    alternativeClusters(fx > cluster1Bounds(end)) = 0;

    
end

