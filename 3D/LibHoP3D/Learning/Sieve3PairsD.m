% this is to filter out rows representing elements with depth
% discontinuities

% cluster3Depths(n2Clusters, n2Clusters, 4, 2)

function [inds, statistics, curTS] = Sieve3PairsD(statistics, curTS, cluster3Depths, lenDisp)
    
    disp('Sieve statistics (step 2)...');
    inds = [];
    
    for i = 1:curTS
        is_ok = true;
        curLine = statistics(i, :);
        center = curLine(1);
        for j = 1:lenDisp
            curDepth =  curLine(2*j+1);
            depthMin =  cluster3Depths(center, curLine(2*j), j, 1);
            depthMax =  cluster3Depths(center, curLine(2*j), j, 2);
            
            if curDepth < depthMin || curDepth > depthMax || depthMin == -999 || depthMax == -999  
                is_ok = false;
                break;
            end
        end
        if is_ok == true
            inds = [inds; i];
        end
    end
    
    statistics = statistics(inds, :);
    curTS = size(statistics, 1);
end
