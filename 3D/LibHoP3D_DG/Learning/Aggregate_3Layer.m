% this is to aggregate statistics for the third layer and prepare it for
% optimization function
% triples should look like that: left, middle, right

function [X, frequencies, triples, is_sparse] = Aggregate_3Layer(statistics, nPrevClusters, curTS, sieve_thresh, is_GPU_USED) 
     
    % [el_cent    el_left depth_left,  el_right depth_right] ->  [el_left  el_cent  el_right]    

    % now we do not care about depths
    nPrevClusters = nPrevClusters + 1;  % to tacle empty cells
    
    
    is_sparse = false;
    if nPrevClusters > 700
        is_sparse = true;
    end
    
    is_GPU_USED = false;
    disp('Warning!!! Gpu is not used in Aggregate_3Layer function');
    
    tic
    disp('create triples matrix');
    % triples: left, center, right
    if is_sparse
        parfor i = 1:nPrevClusters
            triples{i} = sparse(zeros(nPrevClusters, nPrevClusters));
        end
    else  % if is_sparse
        if is_GPU_USED
            triples = gpuArray.zeros(nPrevClusters, nPrevClusters, nPrevClusters);
        else
            triples = zeros(nPrevClusters, nPrevClusters, nPrevClusters);
        end
    end
    

    disp('fill triples matrix');
    
    if ~is_GPU_USED    % everything is done on CPU
        if ~is_sparse
            % left center right
            for i = 1:curTS
                triples(statistics(i,1), statistics(i,2), statistics(i,3)) = triples(statistics(i,1), statistics(i,2), statistics(i,3)) + 1;
            end
        else  % for
            for i = 1:curTS
                triples{statistics(i,1)}(statistics(i,2), statistics(i,3)) = triples{statistics(i,1)}(statistics(i,2), statistics(i,3)) + 1;
            end
        end
        
    elseif is_GPU_USED   % everything is done on GPU

        if ~is_sparse

            ind = sub2ind([nPrevClusters, nPrevClusters, nPrevClusters], int32(statistics(:,1)), int32(statistics(:,2)), int32(statistics(:,3)));
            ind = uint32(ind);
            indG = gpuArray(ind);
            edges = 1:1:(nPrevClusters*nPrevClusters*nPrevClusters);
            hh = histc(indG, edges);

            iinds = hh > 0;
            edges = edges(iinds);
            triples(edges) = hh(iinds);
            triples = gather(triples);

            clear('edges');
        else   % sparse

            for i = 1:curTS
                 triples{statistics(i,1)}(statistics(i,2), statistics(i,3)) = triples{statistics(i,1)}(statistics(i,2), statistics(i,3)) + 1;
            end
        end
    end
    
    
    tic
    disp('aggregate statistics');
    ind = 0;
    X = zeros(curTS,3);
    frequencies = zeros(curTS,1);

    for i = 1:nPrevClusters
        for j = 1:nPrevClusters
            for k = 1:nPrevClusters
                if ~is_sparse
                    if triples(i,j,k) >= sieve_thresh
                        ind = ind + 1;

                        X(ind,1) = i;  % left   or top
                        X(ind,2) = j;  % centre  or centre
                        X(ind,3) = k;  % right   or bottom
                        frequencies(ind) = triples(i,j,k); 
                    end
                else  % if is_sparse
                    if triples{i}(j,k) >= sieve_thresh
                        ind = ind + 1;

                        X(ind,1) = i;  % left   or top
                        X(ind,2) = j;  % centre  or centre
                        X(ind,3) = k;  % right   or bottom
                        frequencies(ind) = triples{i}(j,k); 
                    end
                end
            end
        end
    end

    % initialization
    X = X(1:ind, :);
    frequencies = frequencies(1:ind);
    
    if is_GPU_USED
        reset(gpuDevice(1));
    end

end




%     X = zeros(curTS,6);
%     table = zeros(2,nPrevClusters);
%     for i = 1:nPrevClusters
%         [clusterX, clusterY] = compute2derivatives(i, nClusters);
%         table(1,i) = clusterX;
%         table(2,i) = clusterY;
%     end

%                     X(ind,1) = table(1, i);
%                     X(ind,2) = table(2, i); % left
%                     X(ind,3) = table(1, j);
%                     X(ind,4) = table(2, j); % center
%                     X(ind,5) = table(1, k);
%                     X(ind,6) = table(2, k); % right

