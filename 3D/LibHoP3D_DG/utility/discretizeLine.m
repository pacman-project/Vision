
% this function convert string of filter responces to the string of
% clucters and errors


function [nearestClusters, errors, alternativeClusters, alternativeErrors] = discretizeLine(fx, strLen, nClusters, cluster1Centres, cluster1Lengths, thresh)

    % compute errors for each point. Errors are interprited as distance to
    % the closest cluster centre edvided by size of the cluster. Distances
    % to the first and last clusters (e.g. 1 and 9) are treated more carefully
    
    errors = zeros(1, strLen);
    nearestClusters = zeros(1, strLen);
    alternativeClusters = zeros(1, strLen);
    alternativeErrors = zeros(1, strLen);
    
    for i = 1:strLen  % for every point
        
        
        %dX <= -thresh +  clusterSizes(1)/2 && dX >= -thresh - clusterSizes(1)/2
        
        % first define two nearest clusters
        if fx(i) < -thresh - cluster1Lengths(1)/2 || fx(i) > thresh + cluster1Lengths(end)/2
            nearestClusters(i) = 0; 
            % everithing else is also 0
            
        elseif fx(i) <= -thresh +  cluster1Lengths(1)/2 && fx(i) >= -thresh - cluster1Lengths(1)/2
            errors(i) = abs(  (fx(i) - cluster1Centres(1))/cluster1Lengths(1)  );
            nearestClusters(i) = 1;
            if  fx(i) < -thresh
                alternativeClusters(i) = 0;
            else
                alternativeClusters(i) = 2;
            end
            if errors(i) > 1
                errors(i) = 1;
            end
            alternativeErrors(i) = errors(i);
            
        elseif fx(i) >= thresh - cluster1Lengths(end)/2 && fx(i) <= thresh + cluster1Lengths(end)/2
            errors(i) = abs(  (fx(i) - cluster1Centres(nClusters))/cluster1Lengths(nClusters)  );
            nearestClusters(i) = nClusters;
            if  fx(i) > thresh
                alternativeClusters(i) = 0;
            else
                alternativeClusters(i) = nClusters - 1;
            end
            alternativeClusters(i) = nClusters;
            if errors(i) > 1
                errors(i) = 1;
            end
            alternativeErrors(i) = errors(i);
        else % the value is a weighted sum of two clusters (main case)

            % here we define left and right clusters
            clusterI = define1Cluster(fx(i), nClusters, cluster1Lengths, thresh);
            errors(i) = abs(  (fx(i) - cluster1Centres(clusterI)) / cluster1Lengths(clusterI)  );
            nearestClusters(i) = clusterI;
            
            % here we define an alternative cluster
            if fx(i) == cluster1Centres(clusterI)
                alternativeClusters(i) = clusterI;
                alternativeErrors(i) = errors(i);
            elseif fx(i) < cluster1Centres(clusterI)
                alternativeClusters(i) = clusterI - 1;
                alternativeErrors(i) = abs(  (fx(i) - cluster1Centres(clusterI - 1)) / cluster1Lengths(clusterI-1)  );
            elseif fx(i) > cluster1Centres(clusterI)
                alternativeClusters(i) = clusterI + 1;
                alternativeErrors(i) = abs(  (fx(i) - cluster1Centres(clusterI + 1)) / cluster1Lengths(clusterI+1)  );
            end
        end      
    end 

end

