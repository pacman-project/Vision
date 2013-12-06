
% this function convert string of filter responces to the string of
% clucters and errors


function [nearestClusters, errors] = discretizeLine(fx, strLen, nClusters, cluster1Centres, cluster1Lengths, thresh)

    % compute errors for each point. Errors are interprited as distance to
    % the closest cluster centre edvided by size of the cluster. Distances
    % to th first and last clusters (e.g. 1 and 9) are treated more carefully
    
    errors = zeros(1, strLen);
    nearestClusters = zeros(1, strLen);
    
    for i = 1:strLen  % for every point
        
        % first define two nearest clusters
        if fx(i) <= cluster1Centres(1)
            errors(i) = abs(  2 *(fx(i) - cluster1Centres(1))/cluster1Lengths(1)  );
            nearestClusters(i) = 1;
            if errors(i) > 1
                errors(i) = 1;
            end
        elseif fx(i) >= cluster1Centres(nClusters)
            errors(i) = abs(  2 * (fx(i) - cluster1Centres(nClusters))/cluster1Lengths(nClusters)  );
            nearestClusters(i) = nClusters;
            if errors(i) > 1
                errors(i) = 1;
            end
        else % the value is a weighted sum of two clusters (main case)

            % here we define left and right clusters
            clusterI = define1Cluster(fx(i), nClusters, cluster1Lengths, thresh);
            errors(i) = abs(  (fx(i) - cluster1Centres(clusterI)) / cluster1Lengths(clusterI)  );
            nearestClusters(i) = clusterI;
        end      
    end 

end

