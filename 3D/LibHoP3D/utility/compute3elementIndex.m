% the function computes element given el

% el = [clX_left, clY_left, clX_middle, clY_middle, clX_right, clY_right];
% triples3Out - list of third layer elements (also 6 dimensional)

function [clusterIndex, is_exist] = compute3elementIndex(el, triples3Out, n3Clusters)

    is_exist = false;
    clusterIndex = 0;
    ind = 0;
    
    while ind < n3Clusters
        ind = ind + 1;
        curEl = triples3Out(ind,:);
        if curEl == el % found
            clusterIndex = ind;
            is_exist = true;
            ind = n3Clusters; % to break the loop
        end
    end
    
end