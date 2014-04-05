% this is to aggregate statistics for the third layer and prepare it for
% optimization function
% triples should look like that: left, middle, right

function [X, frequencies, curTS, triples] = Aggregate_3Layer(statistics, nClusters, n2Clusters, curTS, sieve_thresh) 
     
    % [el_cent    el_left depth_left,  el_right depth_right] ->  [el_left  el_cent  el_right]    

    % now we do not care about depths
    inds3 = [2,1,4];
    statistics = statistics(:, inds3);
    triples = uint16(zeros(n2Clusters, n2Clusters, n2Clusters)); % x-direction

    for i = 1:curTS
        curline = statistics(i,:);
        % left center right
        triples(curline(1), curline(2), curline(3))  = triples(curline(1), curline(2), curline(3)) + 1;
        %top center bottoom
        if mod(i,200000) == 0
            i
        end
    end
    
    table = zeros(2,n2Clusters);
    for i = 1:n2Clusters
        [clusterX, clusterY] = compute2derivatives(i, nClusters);
        table(1,i) = clusterX;
        table(2,i) = clusterY;
    end
    
    ind = 0;
    %X = zeros(curTS,6);
    X = zeros(curTS,3);
    frequencies = zeros(curTS,1);

    for i = 1:n2Clusters
        for j = 1:n2Clusters
            for k = 1:n2Clusters
                if triples(i,j,k) > sieve_thresh
                    ind = ind + 1;
%                     X(ind,1) = table(1, i);
%                     X(ind,2) = table(2, i); % left
%                     X(ind,3) = table(1, j);
%                     X(ind,4) = table(2, j); % center
%                     X(ind,5) = table(1, k);
%                     X(ind,6) = table(2, k); % right

                    X(ind,1) = i;  % left
                    X(ind,2) = j;  % centre
                    X(ind,3) = k;  % right
                    frequencies(ind) = triples(i,j,k); 
                end
            end
        end
    end

    % initialization
    X = X(1:ind, :);
    frequencies = frequencies(1:ind);
    curTS = ind;
end

