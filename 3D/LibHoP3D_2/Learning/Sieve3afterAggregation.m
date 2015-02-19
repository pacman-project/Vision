% sieves statistics after aggregation

function [indsN, statistics] = Sieve3afterAggregation(statistics, triples, sieve_thresh, is_sparse, is_GPU_USED)

    disp('Sieve statistics after aggregation ...');
%     curTS = size(statistics, 1);
%     inds = zeros(1, curTS);
%     
%     lenX = size(X, 1);
    
%     parfor i = 1:curTS
%         
%         if ~is_sparse
%             if triples(statistics(i,2), statistics(i,1), statistics(i,4)) >= sieve_thresh
%                 inds(i) = 1;
%             end
%         else
%             if triples{statistics(i,2)}(statistics(i,1), statistics(i,4)) >= sieve_thresh
%                 inds(i) = 1;
%             end
%         end
%        
% %        if mod(i,100000) == 0
% %            i 
% %        end
%        
%     end

%     indsN = inds == 1;
%     statistics = statistics(indsN, :);

    
    if ~is_sparse
        
        statistics = double(statistics); % this is because of function sub2ind
        if  is_GPU_USED
            triplesG = gpuArray(triples);
            ind = sub2ind(size(triples), statistics(:,2), statistics(:,1), statistics(:,4));
            indsN = triplesG(ind) >= sieve_thresh;
            statistics = statistics(indsN, :);
        else
            ind = sub2ind(size(triples), statistics(:,2), statistics(:,1), statistics(:,4));
            indsN = triples(ind) >= sieve_thresh;
            statistics = statistics(indsN, :);
        end
        
        statistics = int16(statistics);
        
    else  % sparse :(  Therefore no GPU
        
        curTS = size(statistics, 1);
        inds = zeros(1, curTS);
        
        parfor i = 1:curTS
            if triples{statistics(i,2)}(statistics(i,1), statistics(i,4)) >= sieve_thresh
                inds(i) = 1;
            end
        end

        indsN = inds == 1;
        statistics = statistics(indsN, :);
        
    end
    
end