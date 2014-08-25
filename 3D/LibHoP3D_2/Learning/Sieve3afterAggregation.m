% sieves statistics after aggregation

function [inds, statistics] = Sieve3afterAggregation(statistics, triples, sieve_thresh)

    disp('Sieve statistics after aggregation ...');
    curTS = size(statistics, 1);
    inds = zeros(1, curTS);
    indsN = 1;
    
%     lenX = size(X, 1);
    
    for i = 1:curTS
       freq = triples(statistics(i,2), statistics(i,1), statistics(i,4));
       
       if freq > sieve_thresh
            inds(indsN) = i;
            indsN = indsN+1;
       end
       
       if mod(i,100000) == 0
           i 
       end
       
    end
    
    inds = inds(1:indsN-1);
    statistics = statistics(inds, :);
end