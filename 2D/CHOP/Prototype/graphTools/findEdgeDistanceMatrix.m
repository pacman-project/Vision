function [ edgeIdMatrix, edgeDistanceMatrix, edgeCoords ] = findEdgeDistanceMatrix( edgeQuantize )
    % Find a round area in which we will look for individual relations.
    edgeIdMatrix = zeros(edgeQuantize, 'int32');
    centerPointSingle = (edgeQuantize+1)/2; 
    yArr = repmat(1:edgeQuantize, edgeQuantize, 1);
    yArrAbs = abs(yArr - centerPointSingle);
    xArr = repmat((1:edgeQuantize)', 1, edgeQuantize);
    xArrAbs = abs(xArr - centerPointSingle);
    distances = sqrt(xArrAbs .^2 + yArrAbs .^2);
    validPointMatrix = distances<=floor(centerPointSingle);
    
    % Identify each point in the area as an individual relation.
    relations = find(validPointMatrix);
    [coordX, coordY] = ind2sub(size(validPointMatrix), relations);
    edgeCoords = [coordX, coordY] - floor(centerPointSingle);
    edgeIdMatrix(relations) = 1:numel(relations);

    % Find the distance matrix with respect to pairs of individual
    % relations.
    points = [xArr(relations), yArr(relations)];
    distMatrix = single(pdist2(points, points, 'euclidean'));
    edgeDistanceMatrix = distMatrix / max(max(distMatrix));
end

