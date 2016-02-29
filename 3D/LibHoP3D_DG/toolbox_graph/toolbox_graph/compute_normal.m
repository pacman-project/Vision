function [normal,normalf] = compute_normal(vertex, face, fileNameNormF)


% compute_normal - compute the normal of a triangulation
%
%   [normal,normalf] = compute_normal(vertex,face);
%
%   normal(i,:) is the normal at vertex i.
%   normalf(j,:) is the normal at face j.
%
%   Copyright (c) 2004 Gabriel Peyré

[vertex,face] = check_face_vertex(vertex,face);

nface = size(face,2);
nvert = size(vertex,2);

% unit normals to the faces
normalf = crossp( vertex(:,face(2,:))-vertex(:,face(1,:)), ...
                  vertex(:,face(3,:))-vertex(:,face(1,:)) );
              

lenNormf = sqrt(sum(normalf.^2,1)); 
if min(lenNormf)< eps
    ids = lenNormf < eps;
    normalf(:, ids) = 1;
    lenNormf = sqrt(sum(normalf.^2,1));
end
normalf = normalf ./ repmat( lenNormf, 3,1 );

%% check normal faces by ray tracing

F1= face(1,:); F2 = face(2,:); F3 = face(3,:);
V1 = vertex(:,F1); V2 = vertex(:,F2); V3 = vertex(:,F3);
FposAll = (V1+V2+V3)/3;  % centre of each face

if nargin == 3
    if exist(fileNameNormF, 'file')
        disp('Warning: Face normals are downloaded from the file!!!');
        aa = load(fileNameNormF);
        normalf = aa.normalf;
    else
        normalf = CheckNormalsFace(vertex, face, FposAll, normalf);  % make sure they are all pointing outside the volume
        save(fileNameNormF, 'normalf');
    end
end

% lenVis = 30;
% idsVis = randperm(nface, lenVis);
% VisualizeTriangulation(face, vertex);
% hold on
% for i = 1:lenVis
%     quiver3(FposAll(1, idsVis(i)), FposAll(2,idsVis(i)), FposAll(3,idsVis(i)), normalf(1,idsVis(i)), normalf(2,idsVis(i)), normalf(3,idsVis(i)));
%     hold on
% end

%% compute normals for each vertex

% unit normal to the vertex
normal = zeros(3,nvert);
for i=1:nface
    f = face(:,i);
    for j=1:3
        normal(:,f(j)) = normal(:,f(j)) + normalf(:,i);
    end
end

LenNorm = sqrt(sum(normal.^2,1));
if min(LenNorm)< eps
    ids = LenNorm < eps;
    normal(:, ids) = 1;
    LenNorm = sqrt(sum(normal.^2,1));
end
normal = normal ./ repmat( LenNorm, 3,1 );

% % normalize
% d = sqrt( sum(normal.^2,1) ); d(d<eps)=1;
% normal = normal ./ repmat( d, 3,1 );








% % enforce that the normal are outward
% v = vertex - repmat(mean(vertex,1), 3,1);
% s = sum( v.*normal, 2 );
% if sum(s>0)<sum(s<0)
%     % flip
%     normal = -normal;
%     normalf = -normalf;
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = crossp(x,y)
% x and y are (m,3) dimensional
z = x;
z(1,:) = x(2,:).*y(3,:) - x(3,:).*y(2,:);
z(2,:) = x(3,:).*y(1,:) - x(1,:).*y(3,:);
z(3,:) = x(1,:).*y(2,:) - x(2,:).*y(1,:);
