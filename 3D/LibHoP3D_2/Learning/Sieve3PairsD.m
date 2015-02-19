% this is to filter out rows representing elements with depth
% discontinuities

% cluster3Depths(n2Clusters, n2Clusters, 4, 2)

function [indSS, statistics, curTS] = Sieve3PairsD(statistics, curTS, cluster3Depths, lenDisp, is_GPU_USED)
    
     disp('Sieve statistics (step 2)...');
     
     if ~is_GPU_USED
        inds = false(1, curTS);

        parfor i = 1:curTS
            is_ok = true;
            curLine = statistics(i, :);

            for j = 1:lenDisp

                curDepth =  curLine(2*j+1);
                depthMin = cluster3Depths(curLine(1), curLine(2*j), j, 1);
                depthMax = cluster3Depths(curLine(1), curLine(2*j), j, 2);

                if curDepth < depthMin || curDepth > depthMax 
                    is_ok = false;
                end
            end

            if is_ok
                inds(i) = 1;
            end
        end

        indSS = inds == 1;
        statistics = statistics(indSS, :);
        curTS = size(statistics, 1);
        
     elseif is_GPU_USED  % compute on GPU

        statisticsG = gpuArray(statistics);
        cluster3DepthsG = gpuArray(cluster3Depths);
        inds = gpuArray.true(curTS, 1);

        for i = 1:lenDisp
            for j = 1:2  % min and max depths

                cluster3DepthsGCur = squeeze(cluster3DepthsG(:,:,i,j));

                indSS  = sub2ind(size(cluster3DepthsGCur), statistics(:, 1), statistics(:, 2*i));
                if j == 1  % min
                    inds = inds & statisticsG(:,2*i+1) >= cluster3DepthsGCur(indSS);
                else  % max
                    inds = inds & statisticsG(:,2*i+1) <= cluster3DepthsGCur(indSS);
                end

            end
        end   

        statisticsG = statisticsG(inds, :);
        statistics = gather(statisticsG);
        indSS = gather(inds);
        curTS = size(statistics, 1);
        
     end
    

end



