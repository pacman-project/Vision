function [absErr, Norm, is_ok] = PlanarTestParaboloid(xs, ys, zs)

    xs2 = xs.^2;
    ys2 = ys.^2;
    xys = xs.*ys;

    X = [xs2', ys2', xys', xs', ys', ones(length(xs),1)];
    if rank(X) < 6
        is_ok = false;
        Norm = NaN;
        absErr = 90;
        return;
    end
    
    is_ok = true;
    A = X\zs';  % (a, b, c, d, e, f)
    
    w = warning('query','last');
    if strcmp(w.identifier, 'MATLAB:singularMatrix');
        a = 2;
    end

    % estimate the normal in the centre of the receptive field (x = 0, y = 0)
    tangentX = [1,0,A(4)];
    tangentY = [0,1,A(5)];
    Norm = cross(tangentX, tangentY);

    Norm = Norm/norm(Norm);
    
    % estimate normal at the point with certain displacement in parameter
    % space
    dists = sqrt(xs2 + ys2 + zs.^2);
    [dists, ids] = sort(dists, 'descend');
    ids = ids(1:round(0.5 * length(dists)));
    
    absErr = 0;
    len = length(ids);
    
    for i = 1:len
        disp = [xs(ids(i)), ys(ids(i))];

        tangentX_disp = [1,0, 2*A(1)*disp(1) + A(3)*disp(2) + A(4)];
        tangentY_disp = [0,1, 2*A(2)*disp(2) + A(3)*disp(1) + A(5)];
        Norm_disp = cross(tangentX_disp, tangentY_disp);
        Norm_disp = Norm_disp/norm(Norm_disp);

        angle = acos(dot(Norm, Norm_disp))*180/pi;
        
        absErr = absErr + abs(angle);
    end
    
    absErr = absErr/len; 
end

