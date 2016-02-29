% root = 'D:\GeodesicMeshTrial\';
% addpath([root,'toolbox_signal']);
% addpath([root,'toolbox_general']);
% addpath([root,'toolbox_graph']);
% addpath([root,'toolbox_wavelet_meshes']);
% addpath([root,'toolbox_fast_marching']);

%  select3dtool

% name = 'D:\GeodesicMeshTrial\toolbox_graph\off\nefertiti.off';
name = 'D:\Input Data\Meshes\Aim@Shape_Selected_1.00\mug20.obj';

[vertex,faces] = read_mesh(name);
clear options; options.name = name;

clf;
options.lighting = 1;
plot_mesh(vertex,faces,options);
shading('faceted');

% nsub = 1; % number of subdivision steps
% options.sub_type = 'loop';
% options.verb = 0;
% [vertex,faces] = perform_mesh_subdivision(vertex0,faces0,nsub,options);



n = size(vertex,2);

clf;
options.lighting = 1;
plot_mesh(vertex,faces,options);
shading('faceted');

dotp = @(u,v)sum(u.*v,1);
R = @(u)reshape(u, [1 1 length(u)]);
Inv1 = @(M,d)[M(2,2,:)./d -M(1,2,:)./d; -M(2,1,:)./d M(1,1,:)./d];
Inv  = @(M)Inv1(M, M(1,1,:).*M(2,2,:) - M(1,2,:).*M(2,1,:));
Mult = @(M,u)[M(1,1,:).*u(1,1,:) + M(1,2,:).*u(2,1,:);  M(2,1,:).*u(1,1,:) + M(2,2,:).*u(2,1,:)];

W = ones(n,1);

I = [500];

U = zeros(n,1);
i = [faces(1,:) faces(2,:) faces(3,:) ];
j = [faces(2,:) faces(3,:) faces(1,:) ];
k = [faces(3,:) faces(1,:) faces(2,:) ];

x  = vertex(:,i);
x1 = vertex(:,j) - x;
x2 = vertex(:,k) - x;

C = [R(dotp(x1,x1)) R(dotp(x1,x2)); R(dotp(x2,x1)) R(dotp(x2,x2))];
S = Inv(C);
a = sum(sum(S));
d11 = sqrt(dotp(x1,x1));
d11 = d11(:);
d22 = sqrt(dotp(x2,x2));
d22 = d22(:);

niter = 3000;
err = zeros(1, niter);

for tt = 1:niter

    uj = U(j);
    uk = U(k);
    u = [R(uj); R(uk)];
    w = R( W(i) );

    b = dotp( sum(S,2), u );
    c = dotp( Mult(S,u), u ) - w.^2;
    delta = max( b.^2 - a.*c, 0);
    d = (b + sqrt(delta) )./a;
    alpha = Mult( S, u - repmat(d, 2, 1) );

    J = find( alpha(1,1,:)>0 | alpha(2,1,:)>0 );
    d1 = d11.*w(:) + uj(:);
    d2 = d22.*w(:) + uk(:);
    d = d(:);
    d(J) = min(d1(J), d2(J));

    U1 = accumarray(i', d, [n 1], @min);  
    U1(U1==0) = Inf;
    U1(I) = 0;
    err(tt) = sum(abs(U(:) - U1(:)));
    U = U1;
    
    if mod(tt, 50) == 0
        disp(tt);
    end
end

if mod(tt, 10) == 0
    figure;
    mycolor = @(U,k)mod(U/max(U), 1/k);
    mycolor = @(U,k)cos(2*pi*k*U/max(U));

    clf;
    options.face_vertex_color = mycolor(U, 5);
    plot_mesh(vertex,faces,options);
    colormap jet(256);
end

k = 10;
pend = round(rand(k,1)*n)+1;
options.method = 'continuous';
paths = compute_geodesic_mesh(U, vertex, faces, pend, options);

clf;
plot_fast_marching_mesh(vertex,faces, mycolor(U, 5), paths, options);












% % name = 'D:\Input Data\Meshes\Aim@Shape_Selected_1.00\mug10.obj'; % 1k
% 
% [V, F] = read_mesh(name);
% VisualizeTriangulation(F, V);
% n = size(V,2);
% 
% %% this is a trial for subdivision of some faces
% 
% % nsub = 1;
% % options.sub_type = 'linear4';
% % options.verb = 0;
% % [V1, F1] = perform_mesh_subdivision(V, F(:, 1:2000), nsub, options);
% % 
% % V = V1;
% % F = [F1, F(:, 2001:end)];
% % figure
% % VisualizeTriangulation(F, V);
% 
% %% this is to measure geodesic distances
% 
% dotp = @(u,v)sum(u.*v,1);
% R = @(u)reshape(u, [1 1 length(u)]);
% 
% Inv1 = @(M,d)[M(2,2,:)./d -M(1,2,:)./d; -M(2,1,:)./d M(1,1,:)./d];
% Inv  = @(M)Inv1(M, M(1,1,:).*M(2,2,:) - M(1,2,:).*M(2,1,:));
% Mult = @(M,u)[M(1,1,:).*u(1,1,:) + M(1,2,:).*u(2,1,:);  M(2,1,:).*u(1,1,:) + M(2,2,:).*u(2,1,:)];
% 
% 
% W = ones(n,1);
