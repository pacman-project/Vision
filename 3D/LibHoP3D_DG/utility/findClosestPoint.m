% finds the closest point from the point to the surface in depth image I

function [x,y,z] = findClosestPoint(I, point)

    [r,c,ch] = size(I);
    
    if point(1) < 0 || point(2) < 0 || point(1) > c || point(2) > r
        disp('Error in findClosestPoint');
    end
    
    [ys, xs] = find(I > 0);
    ind = sub2ind(size(I), ys, xs);
    zs = I(ind);
    
    len = length(zs);
    distX = xs - point(1);
    distY = ys - point(2);
    distZ = zs - point(3);
    
    dists = sqrt(distX.^2 + distY.^2 + distZ.^2);
    idx = find(dists == min(dists));
    
    x = xs(idx(1));
    y = ys(idx(1));
    z = zs(idx(1)); 
end

