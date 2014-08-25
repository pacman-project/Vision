% This is the function to delete cluster
% matirix X (n*d) data samples
% matrix Z (nClusters*d) - cluster centroids
% distances - (n, nCluster) - matrix wuth distances

function  [sampleLables, clusterLables, clusterFrequencies, Z, nClusters] = Isodata_deleteCluster(deletedClusterInds, nDeleted, replaced_marks, sampleLables, clusterLables, clusterFrequencies, Z, nClusters)

    % first make shure deletedClusterInds and are sorted in ascending order
    [deletedClusterInds, ind] = sort(deletedClusterInds, 'ascend');
    replaced_marks = replaced_marks(ind);
    
    % Unmark samples 
    for ii = 1:nDeleted
        curLabel = deletedClusterInds(ii);
        replaceLabel = replaced_marks(ii);
        sampleLables(sampleLables == curLabel) = replaceLabel;
    end

    % shrink matrices Z, clusterFrequencies, clusterLables

    for ii = nDeleted:-1:1
         lineInd = deletedClusterInds(ii);
         Z(lineInd,:) = [];
         clusterFrequencies(lineInd) = [];
         clusterLables(lineInd) = [];
    end

    % reduce nClusters
    nClusters = nClusters - nDeleted;

    % last step is renumeration of marks
    for ii = 1:nClusters
        if clusterLables(ii) ~= ii
            sampleLables(sampleLables == clusterLables(ii)) = ii;
        end
    end  

    clusterLables = 1:1:nClusters;
end