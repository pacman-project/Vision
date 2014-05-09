% sieves statistics after aggregation

function [inds, statistics] = Sieve3afterAggregation(statistics, triples, sieve_thresh, X)

    disp('Sieve statistics after aggregation ...');
    inds = [];
    inds3 = [2,1,4];
    curTS = size(statistics, 1);
    lenX = size(X, 1);
    
    for i = 1:curTS
       curEl = statistics(i,inds3);
       freq = triples(curEl(1), curEl(2), curEl(3));
       
       if freq > sieve_thresh
           inds = [inds, i];
       end
       
       if mod(i,100000) == 0
           i 
       end
       
    end
    
    statistics = statistics(inds, :);
end