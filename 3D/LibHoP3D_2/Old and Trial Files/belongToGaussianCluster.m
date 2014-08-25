% this function shows whether the element belongs to the cluster or not
% it returns the probability

% x has the dimension 3 (the same as mu)

function [p] = belongToGaussianCluster(x, mu, Sigma)
 
    p = sqrt((x-mu)'*inv(Sigma)*(x-mu));

end