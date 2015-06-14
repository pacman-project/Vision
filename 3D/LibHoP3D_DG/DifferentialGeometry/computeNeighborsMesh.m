% this is to compute neighbours from the mesh or point cloud
% V - list of vertices

function  [inds] = computeNeighborsMesh(centrePoint, V, elementRadius)
    
    % bounding box
    indsX = V(:,1) >= centrePoint(1) - elementRadius & V(:,1) <= centrePoint(1) + elementRadius;
    indsY = V(:,2) >= centrePoint(2) - elementRadius & V(:,2) <= centrePoint(2) + elementRadius;
    indsZ = V(:,3) >= centrePoint(3) - elementRadius & V(:,3) <= centrePoint(3) + elementRadius;
    
    inds = indsX & indsY & indsZ;
    
end

