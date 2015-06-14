% returns a list of parts in a receptive field

function [marksOut, indsXOut, indsYOut, depths, dists] = GetPartsNeighbour(I, marks, x, y, z, rad, is2D)

    [r,c,~] = size(I);
    
    if is2D  % compute neigbouring points in just 2D GetPartsNeighbour
        
        disp('Error. Function computeNeighbors is not defined');
        
%         % get a local neighbourhood of the point (circle)
%         [indsXOut, indsYOut] = getDispAbs(2, rad);
%         indsXOut = x + indsXOut;
%         indsYOut = y + indsYOut;
%         
%         % find indices that are out of image
%         inds = indsXOut > 0 & indsXOut <= c & indsYOut > 0 & indsYOut <= r; 
%         indsXOut = indsXOut(inds);
%         indsYOut = indsYOut(inds);
%         
%         ind = sub2ind(size(I), indsYOut, indsXOut);
%         depths = I(ind);
%         inds = find(depths > minDepth & depths < maxDepth);
%         depths = depths(inds);
%         indsXOut = indsXOut(inds);
%         indsYOut = indsYOut(inds);
        
    else % compute neigbouring points in 3D (sphere)
        
        % get a local neighbourhood in 2D (circle)
        [indsXOut, indsYOut] = getDispAbs(2, rad);
        indsXOut = x + indsXOut;
        indsYOut = y + indsYOut;
        
        % find indices that are out of image
        inds = indsXOut > 0 & indsXOut <= c & indsYOut > 0 & indsYOut <= r; %  & dist > 0.8 * rad
        
        [indsXOut, indsYOut] = myFilter2(indsXOut, indsYOut, inds);      
%         indsXOut = indsXOut(inds);
%         indsYOut = indsYOut(inds);
        
        ind = sub2ind(size(I), indsYOut, indsXOut);
        depths = I(ind);
        
        % now check which of these points belong to a sphere centered at [x,y,z]
        xs = indsXOut - x;
        ys = indsYOut - y;
        zs = depths - z;
        
        % compute a distance from the centre of the receptive field
        dists = xs.^2 + ys.^2 + zs.^2;
        inds = find(dists < rad^2);
        
        [indsXOut, indsYOut, depths, dists] = myFilter4(indsXOut, indsYOut, depths, dists, inds);
%         indsXOut = indsXOut(inds);
%         indsYOut = indsYOut(inds);
%         depths = depths(inds);
                
        ind = sub2ind(size(I), indsYOut, indsXOut);
        
        % retrieve parts in the neighbourhood

        marksOut = marks(ind);
        inds = find(marksOut > 0);
        
%         if length(inds) > 1
%             inds = inds(1);     % to do something smarter!!!!!!!!
%         end
        
        marksOut = marksOut(inds);
        [indsXOut, indsYOut, depths, dists] = myFilter4(indsXOut, indsYOut, depths, dists, inds);
%         indsXOut = indsXOut(inds);
%         indsYOut = indsYOut(inds);
%         depths = depths(inds);

    end


end

function [af, bf] = myFilter2(a,b, ids)
    af = a(ids);
    bf = b(ids);
end

function [af, bf, cf] = myFilter3(a,b,c, ids)
    af = a(ids);
    bf = b(ids);
    cf = c(ids);
end

function [af, bf, cf, df] = myFilter4(a,b,c,d, ids)
    af = a(ids);
    bf = b(ids);
    cf = c(ids);
    df = d(ids);
end
