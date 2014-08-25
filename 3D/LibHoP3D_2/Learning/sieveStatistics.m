% this function is to sieve statistics of all layers

function [statistics, outputCoords, cluster3Depths, curTS] = sieveStatistics(statistics, outputCoords, numDisps, thresh3Pair)

        curTS = size(statistics, 1);

        % for spead up reasons we sort statistics by the first column (central element)
        [~, order] = sort(statistics(:,1));
        statistics = statistics(order,:);
        outputCoords = outputCoords(order,:);
        clear('order');

        % compute the most frequent pairs
        [tablePairs, numPairs] = Learn3Pairs(statistics, curTS, n2Clusters, thresh3Pair, numDisps);
        
        % filter out rows with the least frequent pairs. those pairs with less than thresh3Pair occurrances will be filtered out
        [ind, statistics, curTS] = Sieve3Pairs(statistics, curTS, tablePairs, numDisps);
        outputCoords = outputCoords(ind, :);
        
        % now we have to measure depth and eliminate rows with depth discontinuities
        [cluster3Depths] = compute3Depths(statistics, n2Clusters, quant, numDisps);
        
        % now we have to filter out the depths with discontinuity
        [ind, statistics, curTS] = Sieve3PairsD(statistics, curTS, cluster3Depths, numDisps);
        outputCoords = outputCoords(ind,:);

end

