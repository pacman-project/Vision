% this function computes centre of gravity of each triangle
%   Copyright (c) 2015   Vladislav Kramarev

function [Vadd, Nadd] = compute_gravitycentre(V, F, N)

    Vadd = (V(:,F(1,:)) + V(:,F(2,:)) + V(:,F(3,:))) / 3;
    Nadd = (N(:,F(1,:)) + N(:,F(2,:)) + N(:,F(3,:))) / 3;

end


