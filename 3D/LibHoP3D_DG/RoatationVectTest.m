% this is to try rotation of one vector a to b

% http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/897677#897677
% answer 14 from there

GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;  norm(cross(A,B)) dot(A,B)  0;  0 0 1]; % a 2D roatational matrix
FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
UU = @(Fi,G) Fi*G*inv(Fi);
U = UU(FFi(a,b), GG(a,b));
RotMatr = @(a,b) UU(FFi(a,b), GG(a,b));

a = rand(3,1); 
b = rand(3,1);

a = a/norm(a);
b = b/norm(a);
R = RotMatr(a,b);