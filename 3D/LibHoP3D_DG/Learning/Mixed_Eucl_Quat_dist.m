% computes a mixed distance between sets of points X and Z
% distance is a weighted sum of Euclidian and quaternion distances 

% matirix X (n*d) data samples
% matrix Z (nCl*d) - cluster centroids
% distances - (n, nCluster) - matrix wuth distances

function distances = Mixed_Eucl_Quat_dist(X, Z, alpha)

    if size(X, 1) ~= 6
        X = X';
    end
    if size(Z, 1) ~= 6
        Z = Z';
    end

    if nargin == 2
        % define the weigting parameter according to the following
        % proportion: dist(3 degrees) = 0.001
        % that is 0.0185 * alpha = 0.001
        alpha = 0.0541;
    end
    
    distances = pdist2(X(1:3, :)', Z(1:3,:)') + alpha * QuanternionDistanceVectorized(X(4:6,:), Z(4:6,:));
end 
    
        