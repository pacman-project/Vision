% this is to sieve statistics deleting non-frequent pairs

% the function write the results in place, i.e to the variable 
% lenDisp should be equal to 2 in this case

function [inds, statistics, curTS] = Sieve3Pairs(statistics, curTS, tablePairs, lenDisp)
    
    disp('Sieve statistics ...');
    inds = [];
    
    for i = 1:curTS
        is_ok = true;
        curLine = statistics(i,:);
        center = curLine(1);
        for j = 1:lenDisp
            if tablePairs(j, center, curLine(2*j)) == 0;
                is_ok = false;
                break;
            end
        end
        if is_ok
            inds = [inds; i];
        end
    end
    
    statistics = statistics(inds, :);
    curTS = size(statistics, 1);
end