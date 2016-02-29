% computes distance of two quaternions 
% according to ALG 5 of Kuffner's paper

% X and Z are lists of quaternions

function distQ = QuanternionDistanceVectorized1(X, Z)
    
    if size(X, 1) ~= 4
        X = X';
    end
    if size(Z, 1) ~= 4
        Z = Z';
    end
    n = size(X, 2);
    nCl = size(Z,2);
    
    lambda = X' * Z;
    
    distQ = -abs(lambda) + 1;
    distQ(isnan(distQ)) = 1;
end

