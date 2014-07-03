function [errors, alternativeErrors] = computeLineReconstructionErrors(nearestClusters, clustersAlt, cluster1Centres, fx, inds, thresh)

    elementSize = 5;
    middle = 3;
    lenstr = length(fx);
    a = zeros(elementSize, lenstr);
    b = zeros(elementSize, lenstr);
        
    
    % signal values
    a(middle, :) = fx;
    a(1,3:end) = fx(1:end-2);
    a(2,2:end) = fx(1:end-1);
    a(4,1:end-1) = fx(2:end);
    a(5,1:end-2) = fx(3:end);
    
    for i = 1:2  % for the nearestClusters and clustersAlt
        
        if i == 2
            nearestClusters = clustersAlt;
        end
    
        % reconstructed values
        diffLine = zeros(1, lenstr);
        diffLine(inds) = nearestClusters;
        indss = find(diffLine > 0);
        diffLine(indss) = cluster1Centres(diffLine(indss));

        b(middle, :) = fx;
        b(1,:) = fx - 2 * diffLine;
        b(2,:) = fx - 1 * diffLine;
        b(4,:) = fx + 1 * diffLine;
        b(5,:) = fx + 2 * diffLine;

        b(1,1) = 0; b(1,2) = 0; b(2,1) = 0;
        b(4, end) = 0; b(5, end) = 0; b(5, end - 1) = 0;

        % compute error measures

        aShort = a(:, inds);
        bShort = b(:, inds);

        err = (aShort-bShort).^2;
        if i == 1
            errors = sqrt(sum(err, 1)) / thresh;
            errors(nearestClusters == 0) = 0;
        else
            alternativeErrors = sqrt(sum(err, 1)) / thresh;
            alternativeErrors(nearestClusters == 0) = 0;
        end
        
    end

end

