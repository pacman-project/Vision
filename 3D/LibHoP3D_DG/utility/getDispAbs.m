function [indsXOut, indsYOut] = getDispAbs(elementType, maxRadius)

% here we define indexes w.r.t a central point (0,0)
    

    if elementType == 1 % point
        indsX = 0;
        indsY = 0;
    elseif elementType == 2 % disk
        SE = strel('disk', maxRadius, 0);
        a = getneighbors(SE)';
        indsX = a(1,:);
        indsY = a(2,:);
    elseif elementType == 3 % square
        SE = strel('square', maxRadius + 1);
        a = getneighbors(SE)';
        indsX = a(1,:);
        indsY = a(2,:);
    end
    
    % sort inds according to distance to centre
    dists = sqrt(indsX.^2 + indsY.^2); 
    [~, inds] = sort(dists, 'ascend');
    indsXOut = indsX(inds);
    indsYOut = indsY(inds);

end

