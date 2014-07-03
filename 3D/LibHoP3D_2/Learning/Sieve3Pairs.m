% this is to sieve statistics deleting non-frequent pairs

% the function write the results in place, i.e to the variable 
% lenDisp should be equal to 2 in this case

function [inds, statistics, curTS] = Sieve3Pairs(statistics, curTS, tablePairs, lenDisp)
    
    disp('Sieve statistics ...');
    inds = zeros(1, curTS);
    indsN = 1;
    
    for i = 1:curTS
        is_ok = true;
        for j = 1:lenDisp
            if tablePairs(j, statistics(i,1), statistics(i,2*j)) == 0;
                is_ok = false;
            end
        end
        if is_ok
            inds(indsN) = i;
            indsN = indsN+1;
        end
    end
    
    inds = inds(1:indsN-1);
    statistics = statistics(inds, :);
    curTS = size(statistics, 1);
end