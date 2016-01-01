% This is to implement step 5 of Isodata algorithm
% matirix X (n*d) data samples
% matrix Z (nCl*d) - cluster centroids
% distances - (n, nCluster) - matrix wuth distances

function  [delta, deltaOverall] = Isodata_distances_step5_Quaternion(distances, sampleLables, frequencies, clusterFrequencies, n, nCl)

    delta = zeros(nCl,1); % initialization
    denominators = zeros(nCl,1); 
    
    for ii = 1:n
        curCluster = sampleLables(ii);
        if curCluster~=0 % if the pixel is marked
            dist = distances(ii, curCluster);
            delta(curCluster) = delta(curCluster) + dist * frequencies(ii);
            denominators(curCluster) = denominators(curCluster) + frequencies(ii);
        end
    end

    delta = delta./denominators;
    deltaOverall = (sum(clusterFrequencies .* delta))/(sum(clusterFrequencies));
end