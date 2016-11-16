function Y=transformvertices(X,M)
% Transform the vertices with 3x3 or 4x4 matrix M
% 
% Y=transformvertices(X,M);
%
% input,
%    X : The vertices
%    M : The 3x3 or 4x4 transformation matrix
%
% output,
%    Y : The transformed vertices

switch size(M,2)
    case 3
        Y=zeros(size(X));
        Y(:,1)= M(1,1).*X(:,1) + M(1,2).*X(:,2) + M(1,3);
        Y(:,2)= M(2,1).*X(:,1) + M(2,2).*X(:,2) + M(2,3);
    case 4
        Y=zeros(size(X));
        Y(:,1)= M(1,1).*X(:,1) + M(1,2).*X(:,2) + M(1,3).*X(:,3) + M(1,4);
        Y(:,2)= M(2,1).*X(:,1) + M(2,2).*X(:,2) + M(2,3).*X(:,3) + M(2,4);
        Y(:,3)= M(3,1).*X(:,1) + M(3,2).*X(:,2) + M(3,3).*X(:,3) + M(3,4);
end

