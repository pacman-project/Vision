% in this script we test a toolbox toolbox_graph

root = 'F:\LibHoP3D\DifferentialGeometry\';

addpath([root, 'toolbox_graph\toolbox_graph']);
addpath([root, 'toolbox_graph\toolbox_graph\off']);
addpath([root, 'toolbox_graph\toolbox_graph\toolbox']);


% %% mesh loading and displaying
% 
% % load the mesh
% name = 'elephant-50kv';
% options.name = name; % useful for displaying
% [vertex,faces] = read_mesh(name);
% % display the mesh
% clf;
% plot_mesh(vertex, faces);
% shading interp;
% 
% %% compute the normal to the mesh, and displace the position of the vertices along the normal
% options.face_vertex_color =  [];
% vertex1 = perform_normal_displacement(vertex,faces,.03);
% 
% clf;
% subplot(1,2,1);
% plot_mesh(vertex,faces,options); shading interp; axis tight;
% subplot(1,2,2);
% plot_mesh(vertex1,faces,options); shading interp; axis tight;


% %% normals for each vertex and each face
% name = 'mushroom';
% options.name = name; % useful for displaying
% [vertex,faces] = read_mesh(name);
% % compute normal per vertex and per face
% [normal,normalf] = compute_normal(vertex,faces);
% % display
% options.normal = normal;
% clf; plot_mesh(vertex,faces,options); shading interp; axis tight;
% options.normal = [];



%% curvature trial

% name = 'mushroom';
% options.name = name;
% rep = ['results/curvature/' name '/'];
% if not(exist(rep))
%     mkdir(rep);
% end
% 
% 
% [vertex,face] = read_mesh(name);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% local covariance analysis
% options.covariance_smoothing = 15;
% [C,U,D] = compute_mesh_local_covariance(vertex,face,vertex,options);
% 
% % options for display
% tau = 1.2;
% options.normal_scaling = 1.5;
% 
% 
% options.normal = squeeze(U(:,2,:));
% clf;
% options.face_vertex_color = perform_saturation( -D(2,:)' - D(3,:)',tau);
% plot_mesh(vertex,face, options);
% shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-covariance.png'], 'png');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% curvature
% 
% options.curvature_smoothing = 10;
% [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,face,options);
% 
% options.normal = [];
% clf;
% options.face_vertex_color = perform_saturation(Cmax,tau);
% plot_mesh(vertex,face, options);
% shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-cmax.png'], 'png');
% clf;
% options.face_vertex_color = perform_saturation(Cmin,tau);
% plot_mesh(vertex,face, options);
% shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-cmin.png'], 'png');
% clf;
% options.face_vertex_color = perform_saturation(Cmean,tau);
% plot_mesh(vertex,face, options);
% shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-cmean.png'], 'png');
% clf;
% options.face_vertex_color = perform_saturation(Cgauss,tau);
% plot_mesh(vertex,face, options);
% shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-cgauss.png'], 'png');
% clf;
% options.face_vertex_color = perform_saturation(abs(Cmin)+abs(Cmax),tau);
% plot_mesh(vertex,face, options);
% shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-cabs.png'], 'png');


% load the mesh
name = 'elephant-50kv';
options.name = name; % useful for displaying
[vertex,faces] = read_mesh(name);
% compute the curvature
options.curvature_smoothing = 10;
options.verb = 0;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,faces,options);
% display
clf;
subplot(1,2,1);
options.face_vertex_color = perform_saturation(Cgauss,1.2);
plot_mesh(vertex,faces, options); shading interp; colormap jet(256);
title('Gaussian curvature');
subplot(1,2,2);
options.face_vertex_color = perform_saturation(abs(Cmin)+abs(Cmax),1.2);
plot_mesh(vertex,faces, options); shading interp; colormap jet(256);
title('Total curvature');








