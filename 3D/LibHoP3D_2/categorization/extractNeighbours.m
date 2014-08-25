% this is the function to extract a list of nearest elements

function [ el_list, weights, numEl ] = extractNeighbours(el)

    [clusterX, clusterY] = compute2derivatives(el, 9);
    el_list = zeros(1,8);
    weights = zeros(1,8);
    numEl = 0;
    
    if clusterX > 1
        curEl = compute2elementIndex(clusterX - 1, clusterY, 9);
        numEl = numEl + 1;
        el_list(numEl) = curEl;
        weights(numEl) = 2;
    end
        
    if clusterX < 9
        curEl = compute2elementIndex(clusterX + 1, clusterY, 9);
        numEl = numEl + 1;
        el_list(numEl) = curEl;
        weights(numEl) = 2;
    end
    
    if clusterY > 1
        curEl = compute2elementIndex(clusterX, clusterY - 1, 9);
        numEl = numEl + 1;
        el_list(numEl) = curEl;
        weights(numEl) = 2;
    end
    
    if clusterY < 9
        curEl = compute2elementIndex(clusterX, clusterY + 1, 9);
        numEl = numEl + 1;
        el_list(numEl) = curEl;
        weights(numEl) = 2;
    end
    % ---------------------------------------------------------------------
    
    if clusterX > 1 && clusterY > 1
        curEl = compute2elementIndex(clusterX - 1, clusterY - 1, 9);
        numEl = numEl + 1;
        el_list(numEl) = curEl;
        weights(numEl) = 1;
    end
        
    if clusterX < 9 && clusterY < 9
        curEl = compute2elementIndex(clusterX + 1, clusterY + 1, 9);
        numEl = numEl + 1;
        el_list(numEl) = curEl;
        weights(numEl) = 1;
    end
    
    if clusterX > 1 && clusterY < 9
        curEl = compute2elementIndex(clusterX - 1, clusterY + 1, 9);
        numEl = numEl + 1;
        el_list(numEl) = curEl;
        weights(numEl) = 1;
    end
    
    if clusterX < 9 && clusterY > 1
        curEl = compute2elementIndex(clusterX + 1, clusterY - 1, 9);
        numEl = numEl + 1;
        el_list(numEl) = curEl;
        weights(numEl) = 1;
    end
    
    el_list = el_list(1:numEl);
    weights = weights(1:numEl);

end

