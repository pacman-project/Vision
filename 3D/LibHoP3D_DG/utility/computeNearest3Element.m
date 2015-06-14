% the function finds the nearest element of the third layer

% el = [clX_left, clY_left, clX_middle, clY_middle, clX_right, clY_right];
% triples3Out - list of third layer elements (also 6 dimensional)


function [elementIndex, dist] = computeNearest3Element(triples3Out, el, n3Clusters)

    [distances] = Isodata_distances(el, triples3Out, 1, n3Clusters, false, false);
    dist = min(distances);
    indMin = find(distances == dist);
    elementIndex = indMin(1);
end
