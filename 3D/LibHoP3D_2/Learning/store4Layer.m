% this is to store the third layer in a format required for visualization and recognition

% input:
% cluster3Depths = zeros(n2Clusters, n2Clusters, lenDisp, 2); 

% triple3OutDepth has the folowing format:

% input:
% cluster4Depths = zeros(n2Clusters, n2Clusters, lenDisp, 2); 

% triple3OutDepth has the folowing format:
% left3 left_depthMin, left_depthMax, left_avg, central3, right3 right_depthMin, right_depthMax, right_avg

function [triple4OutDepth] = store4Layer(triples4Out, cluster4Depths, n4Clusters, nClusters, partSelectionMethod)

    triple4OutDepth = zeros(n4Clusters, 9); % 4+1+4  % top - central - bottom 
    
    for i = 1:n4Clusters
        cur = triples4Out(i,:);
        
        if partSelectionMethod == 1 % rewrite to format: leftX leftY centralX, centralY, rightX rightY
            cur = [cur(1), cur(2) cur(3)];
        end
        
        line = zeros(1,9);
        
        inds = [1,5,6];
        line(inds) = cur;
        
        % now compute relative depths of the left and right elements
        left = cur(1);
        center = cur(2);
        right = cur(3);
        
        line(2) = cluster4Depths(center, left, 1, 1);
        line(3) = cluster4Depths(center, left, 1, 2);
        line(4) = cluster4Depths(center, left, 1, 3);
        line(7) = cluster4Depths(center, right, 2, 1);
        line(8) = cluster4Depths(center, right, 2, 2);
        line(9) = cluster4Depths(center, right, 2, 3);

        triple4OutDepth(i,:) = line;
    end
    
end