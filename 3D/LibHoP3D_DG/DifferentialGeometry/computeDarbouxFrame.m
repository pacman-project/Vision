% estimate the local frame of reference at this point
% method of Gabriel Taubin

function [V,D] = computeDarbouxFrame(Norm, xs, ys, zs)  % dx, dy,

%     % (1). compute normal of the point

%     tanX = [1, 0, dx];
%     tanY = [0, 1, dy];
%     Norm = cross(tanX, tanY)';
%     Norm = Norm/norm(Norm);


    % (2). estimate a matrix M
    numPoints = length(xs);
    M = zeros(3,3);
    
    % compute weights for each point (reverse to euclidian distances)
    dist = sqrt(xs.^2 + ys.^2 + zs.^2);  
    sd = sum(dist);
    w = dist/sd;

    for ii = 2:numPoints

        vj = [0,0,0]';
        vi = [xs(ii), ys(ii), zs(ii)]';
        vect = (vi - vj);
        if norm(vect) < 1
            continue;
        end

        Tij = (eye(3) - Norm*Norm')*vect;
        Tij = Tij / norm(Tij);

        % estimate directional curvature
        kij = 2 * Norm' * vect / (norm(vect)^2); 

        MTemp = w(ii) * kij * Tij * Tij';

        M = M + MTemp;
    end

    [V,D] = eig(M);
    values = abs(diag(D));
    
    % if this is a planar patch
    if values(1) < 10^-3  && values(2) < 10^-3  && values(3) < 10^-3
        V(:,1) = Norm;
        Xtemp = [1,0,0]';
        V(:,3) = cross(V(:,1), Xtemp);
        V(:,2) = cross(V(:,3), V(:,1));
    else
    
        [valuesSort, idx] = sort(values, 'ascend');
        D = diag(valuesSort);
        V = V(:, idx);

        % if the direction of normal is reverced
        if abs(sum(V(:,1)+ Norm)) < 10^-5
            V(:,1) = -V(:,1);
        end

        % if this is the umbilic point
        if abs(valuesSort(2) - valuesSort(3)) < 10^-3  
            % make principal directions aligned with the main axis
            Xtemp = [1,0,0]';
            V(:,3) = cross(V(:,1), Xtemp);
            V(:,2) = cross(V(:,3), V(:,1));
        end
    end
    
    V(:, 1) = V(:,1)/norm(V(:,1));
    V(:, 2) = V(:,2)/norm(V(:,2));
    V(:, 3) = V(:,3)/norm(V(:,3));

end

