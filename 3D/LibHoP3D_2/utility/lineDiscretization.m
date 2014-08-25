% for this function we assum that centers of clusters are at equal
% distances from each other

function [output] = lineDiscretization(fx, strlen, nClusters, kernelW, clusterCenters, clusterSize, thresh)

    % |----1-----------2----------3-----------4----------5---------6----------7----|
    % or
    % |-----1----------2----------3----------4----------5-----|
    

    filterResp = double(zeros(nClusters, strlen));  % sum of filter responses in the area
    output = zeros(size(fx));
    
    % now we define responces to two nearest filters in every point
    for i = 1:(strlen - 1) % for every pixel
        % first define two nearest clusters
        if fx(i) <= clusterCenters(1)
            filterResp(1, i) = 1.0;
        elseif fx(i) >= clusterCenters(nClusters)
            filterResp(nClusters,i) = 1.0;
        else % the value is a weighted sum of two clusters (main case)

            % here we define left and right clusters
            clusterLeft = define1Cluster(fx(i), nClusters, clusterSize, thresh);
            if fx(i) < clusterCenters(clusterLeft)
                clusterLeft = clusterLeft - 1;
            end
            clusterRight = clusterLeft + 1;

            
            % then we have to compute contribution of both filters
            weightR = 1.0 - (clusterCenters(clusterRight) - fx(i)) / clusterSize(clusterRight);
            weightL = 1.0 - (fx(i) - clusterCenters(clusterLeft))  / clusterSize(clusterLeft);
            % write these weights to the matrix
            filterResp(clusterRight, i) = weightR;
            filterResp(clusterLeft, i) = weightL;
        end      
    end   % it should work perfectly

    % filter responces averaged over the cirtain area
    kernel = [1.0, 1.0, 1.0, 1.0, 1.0]; % kernel = [0.6, 0.8, 1.0, 0.8, 0.6];
    filterRespA = zeros(nClusters, strlen); 
   % coverage = zeros(nClusters, strlen);
    for i = 1:nClusters
        line = filterResp(i,:);
        %linef = conv(line, kernel, 'same');
        linef = conv(line, kernel);
        linef = linef(3:end - 2);  % it should work instead
        filterRespA(i, :) = linef;
    end

    % performing discretization step here
    cur = 0;
    while (cur < strlen - 7)
        cur = cur + 3;
        % look for the max value in three rows (cur, cur+1, cur+2, ... , cur+5)
        matr =  filterRespA(:,cur:cur+4);
        % kernelW = [0.0756, 0.2127, 0.2834, 0.2127, 0.0756];
        for rr = 1:nClusters
            matr(rr,:) = matr(rr,:) .* kernelW;
        end

        [r,c] = find(matr == max(max(matr)));
        % if there are two similar maximums
        if length(r) > 1
            r = r(length(r));
            c = c(length(c));  % for a moment
        end

        col = cur + c - 1;
        % mark = [1,3,5,3,1];
        % coverage(r, col-2:col+2) = mark;
        output(col) = r;
        cur = col; 
        
        % TODO: here we might also do something about the error measure.
    end 
    

    

end