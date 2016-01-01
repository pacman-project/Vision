
function [partIDs, nParts] = symmetricTrial(nClusters)

    nParts = (nClusters + 1) * nClusters/2;
    d = 1:nParts;
    partIDs = zeros(nClusters, nClusters);
    cur = 1;
    
    for i = 1:nClusters
        len = nClusters - i;
        partIDs(i, i:i+len) = d(cur:cur + len);
        partIDs(i:i+len, i) = d(cur:cur + len)';
        cur = cur + len + 1;
    end
end

