% This function draws cluster in 3D
% Given mu and covariance matrix Sigma

function [out] = drawCluster(mu, Sigma)

    scales = [0.5, 1.0, 1.5, 2.0]; % radius in terms of sigma
    len = length(scales);
    
    % L: eigenvalue diagonal matrix
    % U: eigen vector matrix, each column is an eigenvector
    [U,L] = eig(Sigma);
    
    for i = 1:len
        hold on
        scale = scales(i);
        radii = 1.5*scale*sqrt(diag(L));
        [xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));
        
        % rotate data with orientation matrix U and center mu
        a = kron(U(:,1),xc); 
        b = kron(U(:,2),yc); 
        c = kron(U(:,3),zc);
        data = a+b+c; n = size(data,2);
        x = data(1:n,:)+mu(1); 
        y = data(n+1:2*n,:)+mu(2); 
        z = data(2*n+1:end,:)+mu(3);
        
        curAlpha = 1 - scale/2.8;
        
        surf(x, y, z, 'FaceColor','red', 'EdgeColor','none');
        camlight 'headlight'; lighting gouraud
        alpha(curAlpha)
    end
    out = true;
    
end

