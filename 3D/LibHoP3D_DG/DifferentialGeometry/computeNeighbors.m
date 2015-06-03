% compute a local neighborhood of the point in 2D or 3D
% range image I
% coordinates

function [indsXOut, indsYOut, depths] = computeNeighbors(I, x, y, rad, is2D, minDepth, maxDepth)

    [r,c,~] = size(I);
    if is2D  % compute neigbouring points in just 2D  
        
        % get a local neighbourhood of the point (circle)
        [indsXOut, indsYOut] = getDispAbs(2, rad);
        indsXOut = x + indsXOut;
        indsYOut = y + indsYOut;
        
        % find indices that are out of image
        inds = indsXOut > 0 & indsXOut <= c & indsYOut > 0 & indsYOut <= r; 
        indsXOut = indsXOut(inds);
        indsYOut = indsYOut(inds);
        
        ind = sub2ind(size(I), indsYOut, indsXOut);
        depths = I(ind);
        inds = find(depths > minDepth & depths < maxDepth);
        depths = depths(inds);
        indsXOut = indsXOut(inds);
        indsYOut = indsYOut(inds);
        
    else % compute neigbouring points in 3D (sphere)
        
        % get a local neighbourhood in 2D (circle)
        [indsXOut, indsYOut] = getDispAbs(2, rad);
        
%         % trial
%         dist = sqrt(indsXOut.^2 + indsYOut.^2);
        
        indsXOut = x + indsXOut;
        indsYOut = y + indsYOut;
        
        % find indices that are out of image
        inds = indsXOut > 0 & indsXOut <= c & indsYOut > 0 & indsYOut <= r; %  & dist > 0.8 * rad
        indsXOut = indsXOut(inds);
        indsYOut = indsYOut(inds);
        
        ind = sub2ind(size(I), indsYOut, indsXOut);
        depths = I(ind);
        
        % now check which of these points belong to a sphere
        xs = indsXOut - x;
        ys = indsYOut - y;
        zs = depths - I(y,x);
        
        % compute a distance from the centre 
        dists = xs.^2 + ys.^2 + zs.^2;
        inds = find(dists < (rad)^2);
        
        depths = depths(inds);
        indsXOut = indsXOut(inds);
        indsYOut = indsYOut(inds);
        
    end


end

