function [t] = rayTriangleIntersectionVectorized (o, d, p0, p1, p2)
% Ray/triangle intersection using the algorithm proposed by Möller and Trumbore (1997).
%
% Input:
%    o : origin.
%    d : direction.
%    p0, p1, p2: vertices of the triangle.
% Output:
%    flag: (0) Reject, (1) Intersect.
%    u,v: barycentric coordinates.
%    t: distance from the ray origin.
% Author: 
%    Jesus Mena

    epsilon = 0.00001;

    e1 = p1-p0;
    e2 = p2-p0;
    %q  = cross(d,e2);
    q = bsxfun(@cross, d', e2')';

%     a  = e1 * q';  %dot(e1,q); % determinant of the matrix M
    a = sum(e1 .* q, 2);
    ids1 = ~(a>-epsilon & a<epsilon);  %  flag = 0; u= 0; v = 0; t = 0;
    
    a = a(ids1);
    p0 = p0(ids1, :);
    q = q(ids1, :);
    e1 = e1(ids1, :);
    e2 = e2(ids1, :);

    f = a.^(-1);
    s = bsxfun(@plus, o, -p0);   % o -p0 ;
    u = f .* sum(s .* q, 2);   % f*(s*q');
    
    ids2 = u >= 0;
    
    q = q(ids2, :);
    e1 = e1(ids2, :);
    e2 = e2(ids2, :);
    s = s(ids2, :);
    f = f(ids2);
    u = u(ids2);
    
    r = bsxfun(@cross, s, e1);
    v = f.* (r(:, 1) * d(1) + r(:, 2) * d(2) + r(:, 3) * d(3));  %dot(d,r);
    ids3  =  ~(v<0.0 | (u+v)>1.0);
    
    f = f(ids3);
    e2 = e2(ids3, :);
    r = r(ids3, :);

    t = f .* sum(e2 .* r, 2); %dot(e2,r); % verified! 
    t = t(t> epsilon);
    a = 2;
    return;
end

