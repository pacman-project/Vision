% this function is to sieve statistics of all layers

function [statistics, outputCoords, clusterCurDepths, curTS] = sieveStatistics(statistics, outputCoords, numDisps, threshCurPair, nPrevClusters, ...
                                    quant, maxRelDepth, is_GPU_USED)

        curTS = size(statistics, 1);

        if is_GPU_USED
        % sort statistics on CPU
        
          [~, order] = sort(statistics(:,1));
          statistics = statistics(order,:);
          outputCoords = outputCoords(order,:);

        elseif is_GPU_USED
            % sort statistics on GPU 
            
            statistics = gpuArray(statistics);
            outputCoords = gpuArray(outputCoords);

            [~, order] = sort(statistics(:, 1));     
            statistics = statistics(order, :);
            outputCoords = outputCoords(order, :);

            statistics = gather(statistics);
            outputCoords = gather(outputCoords);
        end

        % compute the most frequent pairs
        [tablePairs, numPairs] = Learn3Pairs(statistics, curTS, nPrevClusters + 1, threshCurPair, numDisps);

        % filter out rows with the least frequent pairs. those pairs with less than threshCurPair occurrances will be filtered out
        [ind, statistics, curTS] = Sieve3Pairs(statistics, curTS, tablePairs, numDisps, is_GPU_USED);
        outputCoords = outputCoords(ind, :);
        disp(curTS);

        % now we have to measure depth and eliminate rows with depth discontinuities
        [clusterCurDepths] = compute3Depths(statistics, nPrevClusters + 1, quant, numDisps, maxRelDepth, is_GPU_USED);

        % now we have to filter out the depths with discontinuity
        [ind, statistics, curTS] = Sieve3PairsD(statistics, curTS, clusterCurDepths, numDisps, is_GPU_USED);
        outputCoords = outputCoords(ind, :);

        if is_GPU_USED
            reset(gpuDevice(1));
        end

end

