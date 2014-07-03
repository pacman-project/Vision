% this is to sieve curStatistics according to eligibility of pairs

function  [outputStatistics, curTS] = Sieve4LayerStatistics(outputStatistics, lenF, n2Clusters, quant)
    
    %      9 5 8      
    %      2 1 4      
    %      6 3 7
    
    a = [9,5,8; 2,1,4; 6,3,7];
    
    lenDisp = 2;
    statistics = [];
    for i = 1:3
        cur = a(i,:);   % ex.  9,5,8
        cols = [(cur(2)-1)*2, (cur(1)-1)*2, (cur(1)-1)*2+1, (cur(3)-1)*2, (cur(3)-1)*2+1];
        cols(cols == 0) = 1; % if i == 2 (middle element)
        curStatistics = outputStatistics(:, cols);
        if i ~= 2       % we have to recompute relative depths
            absDepthsCol = (cur(2)-1)*2 + 1;
            absDepths = outputStatistics(:, absDepthsCol);
            curStatistics(:,3) = curStatistics(:,3) - absDepths;
            curStatistics(:,5) = curStatistics(:,5) - absDepths;
        end
        statistics = [statistics; curStatistics];
    end
    
    curTS = size(statistics, 1);
    % compute the most frequent pairs
    thresh3Pair = 0.01 * lenF;
    [tablePairs, numPairs] = Learn3Pairs(statistics, curTS, n2Clusters, thresh3Pair, lenDisp);
    
    % now we have to measure depth and eliminate rows with depth discontinuities
    [cluster3Depths] = compute3Depths(statistics, n2Clusters, quant, lenDisp);
    
    
    for i = 1:3
        cur = a(i,:);   % ex.  9,5,8
        cols = [(cur(2)-1)*2, (cur(1)-1)*2, (cur(1)-1)*2+1, (cur(3)-1)*2, (cur(3)-1)*2+1];
        cols(cols == 0) = 1; % if i == 2 (middle element)
        curStatistics = outputStatistics(:, cols);
        if i ~= 2       % we have to recompute relative depths
            absDepthsCol = (cur(2)-1)*2 + 1;
            absDepths = outputStatistics(:, absDepthsCol);
            curStatistics(:,3) = curStatistics(:,3) - absDepths;
            curStatistics(:,5) = curStatistics(:,5) - absDepths;
        end
        curTS = size(curStatistics, 1);

        % filter out rows with the least frequent pairs. those pairs with less than thresh3Pair occurrances will be filtered out
        [inds, curStatistics, ~] = Sieve3Pairs(curStatistics, curTS, tablePairs, lenDisp);
        outputStatistics = outputStatistics(inds, :);
        curTS = size(outputStatistics, 1);
        % now we have to filter out the depths with discontinuity
        [inds, ~, ~] = Sieve3PairsD(curStatistics, curTS, cluster3Depths, lenDisp);
        outputStatistics = outputStatistics(inds, :);
        curTS = size(outputStatistics, 1);
    end
end