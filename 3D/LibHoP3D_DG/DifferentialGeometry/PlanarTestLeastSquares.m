% this is another function to test if the surface is planar
% z = ax + by + c


function [absErr, Norm, is_ok] = PlanarTestLeastSquares(xs, ys, zs)

    X = [xs', ys', ones(length(xs),1)];
    
    if rank(X) < 3
        is_ok = false;
        Norm = NaN;
        absErr = 90;
        return;
    end
    
    is_ok = true;
    A = X\zs';   % A = (a, b, c)
    
    % to get a normal we have to transfer it to the form ax+by+cz+d = 0
    % then the norm becomes (a,b,c)

    Norm = [-A(1), -A(2), 1];
    Norm = Norm/norm(Norm);
    
    % estimate the error
    
    % reconstruction
    b = X * A;
    b = abs(zs' - b);
    absErr = sum(b)/length(xs); 
end

