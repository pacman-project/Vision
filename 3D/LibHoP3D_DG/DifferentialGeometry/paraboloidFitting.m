% this function uses the local fitting method (second order paraboloid)
% computes canonical forms QF1, QF2, gaussian and mean curvatures, 



function [H, K, QF1, QF2, S, V, D, is_ok, zsLocal] = paraboloidFitting(xs, ys, zs)

    xs2 = xs.^2;
    ys2 = ys.^2;
    xys = xs.*ys;

    X = [xs2', ys2', xys', xs', ys', ones(length(xs),1)];
    
    if rank(X) < 6
        H = 0;
        K = 0;
        QF1 = 0;
        QF2 = 0;
        S = 0;
        V = 0;
        D = 0;
        is_ok = false;
        return;
    end
    A = X\zs';  % (a, b, c, d, e, f)
    
    zsLocal = A(1)*xs2 + A(2)*ys2 +A(3)*xys + A(4)*xs + A(5)*ys + A(6);

    % estimate the normal in the centre of the receptive field (x = 0, y = 0)
    tangentX = [1,0,A(4)];
    tangentY = [0,1,A(5)];
    Norm = cross(tangentX, tangentY);

    Norm = Norm/norm(Norm);
    
    % estimate normal at the point with certain displacement in parameter
    % space
%     disp = [2,0];
%     
%     tangentX_disp = [1,0, 2*A(1)*disp(1) + A(3)*disp(2) + A(4)];
%     tangentY_disp = [0,1, 2*A(2)*disp(2) + A(3)*disp(1) + A(5)];
%     Norm_disp = cross(tangentX_disp, tangentY_disp);
%     Norm_disp = Norm_disp/norm(Norm_disp);
%     
%     angle = acos(dot(Norm, Norm_disp))*180/pi;
%     

    % estimate the first and second fundamental forms
    % we have to compute cos(N,z)
    z = [0,0,1];
    cosNZ = dot(Norm, z);
    QF1 = [1+A(4)^2, A(4)*A(5); A(4)*A(5), 1+A(5)^2];
    QF2 = [2*A(1),2*A(3); 2*A(3), 2*A(2)] * cosNZ;

    E = QF1(1,1); F = QF1(2,1); G = QF1(2,2);
    L = QF2(1,1); M = QF2(2,1); N = QF2(2,2);

    % gaussian and mean curvatures   
    K = (L*N - M^2)/(E*G - F^2);
    H = (E*N + G*L -2*F*M) / (2*(E*G - F^2));
    
    % to find principal curvatures we have to solve the equation
    % k^2 + 2*H*k + K = 0

    a = 1; b = - 2*H; c = K;
    Det = b^2 - 4*a*c;
    k1 = (-b - sqrt(Det))/(2*a);
    k2 = (-b + sqrt(Det))/(2*a);

    % to find principal directions we have to solve equation
    % (EM - FL)t^2 + (EN-GL)t + (FN - GM) = 0

    c = E*M - F*L;
    b = E*N - G*L;
    a = F*N - G*M;

    Det = b^2 - 4*a*c;
    t1 = (-b - sqrt(Det))/(2*a);
    t2 = (-b + sqrt(Det))/(2*a);

    z1 = -(Norm(1) + t1 * Norm(2))/Norm(3);
    z2 = -(Norm(1) + t2 * Norm(2))/Norm(3);

    V1 = [1, t1, z1];
    V2 = [1, t2, z2];
    
    V1 = V1/norm(V1);
    V2 = V2/norm(V2);
    
    V = [V1',V2', Norm'];
    D = eye(3);
    D(1,1) = k1;
    D(2,2) = k2;

%     if abs(K) < zeroThresh
%         K = 0;
%     end
%     if abs(H) < zeroThresh
%         H = 0;
%     end

%     % shape operator
    S = inv(QF1)*QF2;
%     [V2,D2] = eig(S);
%     
%     % now form 3 dimensional matrices 
%     V = zeros(3,3);
%     D = zeros(3,3);
%     
%     V(1:2,2:3) = V2;
%     V(:,1) = Norm;
%     V(3,2) = -(Norm(1)*V2(1,1) + Norm(2)*V2(2,1))/Norm(3);
%     V(3,3) = -(Norm(1)*V2(1,2) + Norm(2)*V2(2,2))/Norm(3);
%     D(2,2) = D2(1,1);
%     D(3,3) = D2(2,2);
%     
%     % Normalize eigenvectors
%     V(:,2) = V(:,2)/ norm(V(:,2));
%     V(:,3) = V(:,3)/ norm(V(:,3));

    is_ok = true;
end
















