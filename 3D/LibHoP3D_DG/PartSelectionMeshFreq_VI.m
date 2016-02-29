% this function performs MDL-based part selection on meshes

% description of the datastructures used here

function [partsOut, lenOut, ORTable] = PartSelectionMeshFreq_VI(X, frequencies, D, layerID)
                      
    % sort according to coverage
    numCand = length(frequencies);
    ORTable = zeros(1, numCand);
    
    freqThres{3} = 100;   % this is something like lenSelected in the previous version
    freqThres{4} = 100;
    freqThres{5} = 100;
    freqThres{6} = 100;
    freqThres{7} = 50;
    freqThres{8} = 50;
    
    quantLI{3} = 0.25;
    quantLI{4} = 0.22;
    quantLI{5} = 0.20;
    quantLI{6} = 0.35;
    quantLI{7} = 0.25;
    quantLI{8} = 0.20;
    
    dd = D(:);
    dd = dd(dd > 0);
    quantThresh = quantile(dd, quantLI{layerID});
    
    lenOut = 0;
    
% %     quantLIs = measureClusterSeparability(D);  % IMPORTANT UPDATE (TO DO!!!)
    
    [frequencies, ids] = sort(frequencies, 'descend');
    X = X(ids, :);    
    D = D(ids, ids);
    partsOut = zeros(size(X));
    
    while frequencies(lenOut+1) > freqThres{layerID}
        lenOut = lenOut + 1;
        partsOut(lenOut, :) = X(lenOut, :);
        ORTable(lenOut) = lenOut;
        
        % get parts with similar geometric properties
        
        % approach 1: just pick parts which are similar according to the threshold value
        idsOR = find(D(lenOut, :) <= quantThresh);

        % approach 2: maximize separability: TODO
%         dists = D(lenOut, :);
%         idsOR = find(D(dists) <= quantThresh);
        
         
        idsOR = idsOR(ORTable(idsOR) == 0);
        
        
        ORTable(idsOR) = lenOut;
        frequencies(idsOR) = 0;
        partsOut(idsOR, :) = X(idsOR, :);
        
        [frequencies, ids] = sort(frequencies, 'descend');
        X = X(ids, :);
        D = D(ids, ids);
        partsOut = partsOut(ids, :);
        ORTable = ORTable(ids);
    end


    % exclude empty spaces from the partsOut
    ids = ORTable ~= 0;
    lenOut = nnz(ids);
    ORTable = ORTable(ids);
    partsOut = partsOut(ids, :);
end






