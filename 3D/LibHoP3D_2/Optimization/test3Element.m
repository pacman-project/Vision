% this is to check if elements are valid or not during traversal

function [el2Ind, is_ok] = test3Element(el, nClusters)

    if min(el) < 1 || max(el) > nClusters
        is_ok = false;
        el2Ind = ones(1,3);
    else
        el2Ind = zeros(1,3);
        for i = 1:3
            clusterX = el(2*i-1);
            clusterY = el(2*i);
            el2Ind(i) = compute2elementIndex(clusterX, clusterY, nClusters);
        end
        is_ok = true;
    end

end

