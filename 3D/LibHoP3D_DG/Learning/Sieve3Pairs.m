% this is to sieve statistics deleting non-frequent pairs

% the function write the results in place, i.e to the variable 
% lenDisp should be equal to 2 in this case

function [inds, statistics, curTS] = Sieve3Pairs(statistics, curTS, tablePairs, lenDisp, is_GPU_USED)
    

     disp('Sieve statistics ...');
     
     if ~is_GPU_USED
        % do it on CPU
        inds = false(1, curTS);
        parfor i = 1:curTS
            is_ok = true;
            for j = 1:lenDisp
                if tablePairs(j, statistics(i,1), statistics(i,2*j)) == 0;
                    is_ok = false;
                    break;
                end
            end
            if is_ok
                inds(i) = 1;
            end
        end
        statistics = statistics(inds, :);
        curTS = size(statistics, 1);
        
        
     elseif is_GPU_USED
        % do it on GPU
        statisticsG = gpuArray(statistics);
        tablePairsG = gpuArray(tablePairs);
        inds = gpuArray.true(curTS, 1);

        for i = 1:lenDisp

            tablePairsGCur = squeeze(tablePairsG(i, :, :));
            indSS  = sub2ind(size(tablePairsGCur), statistics(:, 1), statistics(:, 2*i));
            inds = inds & tablePairsGCur(indSS);

        end   

        statisticsG = statisticsG(inds, :);
        statistics = gather(statisticsG);
        inds = gather(inds);
        curTS = size(statistics, 1); 
     end
    
end