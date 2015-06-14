% this function is to sieve statistics of all layers

function [statistics, outputCoords, outputScales, outputFrames, clusterCurDepths, curTS] = sieveStatistics(statistics, outputCoords, outputScales, outputFrames, numDisps,...
                                                        threshCurPair, nPrevClusters, quant, maxRelDepth, is_GPU_USED)

        curTS = size(statistics, 1);

        if ~is_GPU_USED
        % sort statistics on CPU
        
        [~, order] = sort(statistics(:,1));
        statistics = statistics(order,:);
        outputCoords = outputCoords(order,:);
        outputScales = outputScales(order, :);
        outputFrames = outputFrames(order, :);

        elseif is_GPU_USED
            % sort statistics on GPU 
            
            statistics = gpuArray(statistics);
            outputCoords = gpuArray(outputCoords);
            outputScales = gpuArray(outputScales);
            outputFrames = gpuArray(outputFrames);

            [~, order] = sort(statistics(:, 1));     
            statistics = statistics(order, :);
            outputCoords = outputCoords(order, :);
            outputScales = outputScales(order, :);
            outputFrames = outputFrames(order, :);

            statistics = gather(statistics);
            outputCoords = gather(outputCoords);
            outputScales = gather(outputScales);
            outputFrames = gather(outputFrames);
        end

        % compute the most frequent pairs
        [tablePairs, numPairs] = Learn3Pairs(statistics, curTS, nPrevClusters + 1, threshCurPair, numDisps);

        % filter out rows with the least frequent pairs. those pairs with less than threshCurPair occurrances will be filtered out
        [order, statistics, curTS] = Sieve3Pairs(statistics, curTS, tablePairs, numDisps, is_GPU_USED);
        outputCoords = outputCoords(order, :);
        outputScales = outputScales(order, :);
        outputFrames = outputFrames(order, :);
        disp(curTS);

        % now we have to measure depth and eliminate rows with depth discontinuities
        [clusterCurDepths] = compute3Depths(statistics, nPrevClusters + 1, quant, numDisps, maxRelDepth, is_GPU_USED);

        % now we have to filter out the depths with discontinuity
        [order, statistics, curTS] = Sieve3PairsD(statistics, curTS, clusterCurDepths, numDisps, is_GPU_USED);
        outputCoords = outputCoords(order, :);
        outputScales = outputScales(order, :);
        outputFrames = outputFrames(order, :);

        if is_GPU_USED
            reset(gpuDevice(1));
        end

end

