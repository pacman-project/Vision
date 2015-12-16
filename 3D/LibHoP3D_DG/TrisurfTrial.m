
V = [0,0,0; 1,1,1; 0,0,1; 1,1,0];
F = [1,2,3; 1,2,4];
alpha = 0.8;
% color = [0.3,0.3,0.3,; 0.9,0.9, 0.9];
load('settings/colormapCopper.mat');  % cmap
color = [cmap(50,:); cmap(60, :)];
trisurf(F,V(:,1),V(:,2),V(:,3), 'FaceColor', 'flat', 'FaceVertexCData', color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
