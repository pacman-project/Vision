% This is to implement step 5 of Isodata algorithm
% matirix X (n*d) data samples
% matrix Z (nClusters*d) - cluster centroids
% distances - (n, nCluster) - matrix wuth distances

function  [delta, deltaOverall] = Isodata_distances_step5(X, sampleLables, frequencies, clusterFrequencies, Z, n, nClusters)

    delta = zeros(nClusters,1); % initialization
    denominators = zeros(nClusters,1);
    
    for ii = 1:n
        curCluster = sampleLables(ii);
        if curCluster~=0 % if the pixel is marked
            x1 = X(ii,:);
            x2 = Z(curCluster, :);
            dist = sqrt(sum((x1-x2).^2));
            delta(curCluster) = delta(curCluster) + dist * frequencies(ii);
            denominators(curCluster) = denominators(curCluster) + frequencies(ii);
        end
    end

    delta = delta./denominators;
    deltaOverall = (sum(clusterFrequencies .* delta))/(sum(clusterFrequencies));
end