% this script is the first trial to learn something from the mesh

root = 'F:\LibHoP3D\DifferentialGeometry\';

addpath([root, 'toolbox_graph\toolbox_graph']);
addpath([root, 'toolbox_graph\toolbox_graph\off']);
addpath([root, 'toolbox_graph\toolbox_graph\toolbox']);



filename = 'F:\LibHoP3D\DifferentialGeometry\Models\D00002.obj';

[V,N,F] = simpleObjReader(filename);

trisurf(F,V(:,1),V(:,2),V(:,3), 'FaceColor', [0.3, 0.3, 0.3], 'EdgeColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.4);
light('Position',[-1.0,-1.0,100.0],'Style','infinite');
axis equal;
lighting phong;
hold on

% determine object size
 XMin = min(V(:, 1));
 XMax = max(V(:, 1));
 YMin = min(V(:, 2));
 YMax = max(V(:, 2));
 ZMin = min(V(:, 3));
 ZMax = max(V(:, 3));
 sizeX = XMax - XMin;
 sizeY = YMax - YMin;
 sizeZ = ZMax - ZMin;
 objSize = [sizeX, sizeY, sizeZ];
 
 
 
 % Let's try to start in each point of the mesh
 
 % determine local geometric properties in the neighborhood
 
 nPoints = size(V,1);
 
 %% vocabulary learning procedure
 

     
% 1. deternine a local geometric properties in the neighborhood
% a. via graph toolbox
    
options.curvature_smoothing = 0.10;
options.averaging_type = 'distance';
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(V,F,options);

sum(Cmin)
sum(Cmax)
sum(Cmean)
sum(Cgauss)

indsV = randperm(nPoints);
subsetVis = 20;
indsV = indsV(1:subsetVis);

elementRadius = receptiveField{2};

vectColors = eye(3);
vecLen = 0.1;

for i = 1:subsetVis  % nPoints
   
    curPoint = V(indsV(i), :);
    
    x = curPoint(1);
    y = curPoint(2);
    z = curPoint(3);
    
    % compute local darboux frames everywhere

    VV = [Umin(:, indsV(i)), Umax(:, indsV(i)), Normal(:, indsV(i))];
    
     % plot the eigenvectors in 3D
    for ii = 1:3
        curVect = VV(:, ii);
        curColor = vectColors(ii, :);
        XX = [x, x + vecLen * curVect(1)];
        YY = [y, y + vecLen * curVect(2)];
        ZZ = [z, z + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor);
        hold on
    end
end
 
 a = 2;
 
 clear all;


