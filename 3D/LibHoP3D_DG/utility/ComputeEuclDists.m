% this is to compute the Euclidian distances from centralPoint to points

function [distances] = ComputeEuclDists(centralPoint, points, numDim)

    if nargin == 2
        numDim = 3;
    end
    
    centralPoint = CheckDimension(centralPoint, numDim);
    points = CheckDimension(points, numDim);
    
    distances = sqrt(sum((points - repmat(centralPoint, size(points,1), 1)).^2, 2));

end

function V = CheckDimension(V, correctSize)
    if size(V,2)~=correctSize
        V = V';
    end
    if size(V,2)~=correctSize
        error('An array does not have the correct format!');
    end
end

