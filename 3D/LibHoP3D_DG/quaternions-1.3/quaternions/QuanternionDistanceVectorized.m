% computes distance of two quaternions 
% according to ALG 5 of Kuffner's paper

function distQ = QuanternionDistanceVectorized(X, Z)
    
    if size(X, 1) ~= 3
        X = X';
    end
    if size(Z, 1) ~= 3
        Z = Z';
    end
    n = size(X, 2);
    nCl = size(Z,2);
    
    %% compute quaternions
%     function Q = computeQuaternion(v1, v2)
%     a = cross(v1, v2);
%     w = norm(v1) * norm(v2) + v1'*v2;
%     Q = qnorm([a; w]);
%     end
    NormX = sqrt(sum(abs(X).^2,1));
    NormZ = sqrt(sum(abs(Z).^2,1));
    X(1,:) = X(1,:) ./ NormX; X(2,:) = X(2,:) ./ NormX; X(3,:) = X(3,:) ./ NormX;
    Z(1,:) = Z(1,:) ./ NormZ; Z(2,:) = Z(2,:) ./ NormZ; Z(3,:) = Z(3,:) ./ NormZ;
    X1 = X(1, :);     X2 = X(2, :);     X3 = X(3, :);

    Q1 = zeros(n, nCl); Q2 = zeros(n, nCl); Q3 = zeros(n, nCl); 
    Q4 = X'*Z + 1;
 
    for i = 1:nCl     % cross product
        Q1(:,i) = (X2*Z(3,i) - X3*Z(2,i))';
        Q2(:,i) = (X3*Z(1,i) - X1*Z(3,i))';
        Q3(:,i) = (X1*Z(2,i) - X2*Z(1,i))';
    end
    %% normalization of quaternions
    QNorm = sqrt(Q1.^2 + Q2.^2 + Q3.^2 + Q4.^2);
%     Q1 = Q1./QNorm; 
%     Q2 = Q2./QNorm;
%     Q3 = Q3./QNorm; 
    Q4 = Q4./QNorm;
    
    %% compute distance between quaternions    
    distQ = sqrt(1 - abs(Q4));
end

