% this is to store the third layer in a format required for visualization and recognition

% input:
% cluster3Depths = zeros(n2Clusters, n2Clusters, lenDisp, 2); 

% triple3OutDepth has the folowing format:
% leftX leftY left_depthMin, left_depthMax, left_avg, centralX, centralY, rightX rightY right_depthMin, right_depthMax, right_avg

function [triple3OutDepth] = store3Layer(triples3Out, cluster3Depths, n3Clusters, nClusters)

    triple3OutDepth = zeros(n3Clusters, 12); % 5+2+5
    
    for i = 1:n3Clusters
        cur = triples3Out(i,:);

        [clusterXL, clusterYL] = compute2derivatives(cur(1), nClusters);
        [clusterXC, clusterYC] = compute2derivatives(cur(2), nClusters);
        [clusterXR, clusterYR] = compute2derivatives(cur(3), nClusters);

        cur = [clusterXL, clusterYL, clusterXC, clusterYC, clusterXR, clusterYR];
        line = zeros(1,12);
        
        inds = [1,2,6,7,8,9];
        line(inds) = cur;
        
        % now compute relative depths of the left and right elements
        left = compute2elementIndex(cur(1), cur(2), nClusters);
        center = compute2elementIndex(cur(3), cur(4), nClusters);
        right = compute2elementIndex(cur(5), cur(6), nClusters);
        
        line(3) = cluster3Depths(center, left, 1, 1);
        line(4) = cluster3Depths(center, left, 1, 2);
        line(5) = cluster3Depths(center, left, 1, 3);
        line(10) = cluster3Depths(center, right, 2, 1);
        line(11) = cluster3Depths(center, right, 2, 2);
        line(12) = cluster3Depths(center, right, 2, 3);

        triple3OutDepth(i,:) = line;
    end
    
end