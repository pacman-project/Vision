% this is to filter out rows representing elements with depth
% discontinuities

% cluster3Depths(n2Clusters, n2Clusters, 4, 2)

function [inds, statistics, curTS] = Sieve3PairsD(statistics, curTS, cluster3Depths, lenDisp)
    
    disp('Sieve statistics (step 2)...');
    inds = zeros(1, curTS);
    indsN = 1;
    
    for i = 1:curTS
        is_ok = true;
        for j = 1:lenDisp
            curDepth =  statistics(i, 2*j+1);
            depthMin =  cluster3Depths(statistics(i,1), statistics(i,2*j), j, 1);
            depthMax =  cluster3Depths(statistics(i,1), statistics(i,2*j), j, 2);
            
            if curDepth < depthMin || curDepth > depthMax || depthMin == -999 || depthMax == -999  
                is_ok = false;
                break;
            end
        end
        if is_ok == true
            inds(indsN) = i;
            indsN = indsN+1;
        end
    end
    
    inds = inds(1:indsN-1);
    statistics = statistics(inds, :);
    curTS = size(statistics, 1);
end
