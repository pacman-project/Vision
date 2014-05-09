% this is to aggregate statistics for the third layer and prepare it for
% optimization function
% triples should look like that: left, middle, right

function [X, frequencies, curTS, triples4] = Aggregate_4Layer(statistics, n3Clusters, triples3, curTS, sieve_thresh) 

    % this is represented in 18-dimensional space
     
    % [el_cent    el_top depth_top,  el_bottom depth_bottom] ->  [el_top  el_cent  el_bottom]    

    % now we do not care about depths
    cols = [1,2,4];
    statistics = statistics(:, cols);
    triples4 = zeros(n3Clusters, n3Clusters, n3Clusters); % y-direction

    for i = 1:curTS
        curline = statistics(i,:);
        % left center right
        triples4(curline(2), curline(1), curline(3)) = triples4(curline(2), curline(1), curline(3)) + 1;
        %top center bottoom
        if mod(i,10000) == 0
            i
        end
    end
    
    ind = 0;
    X = zeros(curTS, 18); % at this layer X is 18 dimensional
    frequencies = zeros(curTS,1);

    for i = 1:n3Clusters
        for j = 1:n3Clusters
            for k = 1:n3Clusters
                if triples4(i,j,k) > sieve_thresh
                    ind = ind + 1;
                    line = [triples3(i,:), triples3(j,:), triples3(k,:)]; % 18 dimensions totally
                    X(ind,:) = line;
                    frequencies(ind) = triples4(i,j,k); 
                end
            end
        end
    end

    % initialization
    X = X(1:ind, :);
    frequencies = frequencies(1:ind);
    curTS = ind;
end

