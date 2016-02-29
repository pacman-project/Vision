% mu = [1 1 1 1];
% sigma = [1,0,0,0;  0,1,0,0;  0,0,1,0;  0,0,0,1];

mu = [1 1 1];
sigma = [1,0,0;  0,1,0;  0,0,1];

numPoints = 10000;
X = mvnrnd(mu, sigma, numPoints);
invSigma = inv(sigma);
dists = zeros(1, numPoints);

% measure mahalanobis distances
for i = 1 : numPoints
    
   dists(i) = sqrt((X(i, :) - mu) * invSigma * (X(i, :) - mu)');
   
end

dists = sort(dists, 'descend');