% Here we try the artificially created mesh

root = 'D:\LibHoP3D_DG\';

addpath([root, 'toolbox_graph\toolbox_graph']);
addpath([root, 'toolbox_graph\toolbox_graph\off']);
addpath([root, 'toolbox_graph\toolbox_graph\toolbox']);


elementRadius = 10;
is2D = false;

imSize = 300;
I = zeros(imSize, imSize);

methodId = 1; 

%% here we create a primitive
primType = 5; % cylinder
zeroThresh = 10^(-3);
slant = pi/6;
elevation = 1;

if primType == 1
    r = 70;
    h = 200;
    
    cylCentre = [150, 150, 0];
    cylinderStartX = cylCentre(1) - r;
    cylinderEndX = cylCentre(1) + r;
    cylinderStartY = cylCentre(2) - h/2;
    cylinderEndY = cylCentre(2) + h/2;
    
    for y = cylinderStartY:cylinderEndY
        for x = cylinderStartX:cylinderEndX
            
            xCyl = x - cylCentre(1);
            z = sqrt(r^2 - xCyl^2);
            
            I(y,x) = z;
            
        end
    end
    
elseif primType == 2 % sphere
    r = 100;
    centre = [150, 150, 0];
    startX = centre(1) - r;
    endX = centre(1) + r;
    startY = centre(2) - r;
    endY = centre(2) + r;
    
    for y = startY:endY
        for x = startX:endX
            
            % check if the point belongs to the sphere
            if (x - centre(1))^2 + (y - centre(2))^2 < r^2
                z = sqrt(r^2 - (x - centre(1))^2 - (y - centre(2))^2) + centre(3);
                I(y,x) = z;
            end
        end
    end
elseif primType == 5 % trancated cone
    
    rt = 30;
    rb = 100;
    h =  200;
    
    cylCentre = [150, 150, 0];
    cylinderStartX = cylCentre(1) - rb;
    cylinderEndX = cylCentre(1) + rb;
    cylinderStartY = cylCentre(2) - h/2;
    cylinderEndY = cylCentre(2) + h/2;
    
    for y = cylinderStartY:cylinderEndY
        for x = cylinderStartX:cylinderEndX
            
            rCur = rt + (rb - rt) * (y - cylinderStartY) / h;
            
            if abs(cylCentre(1) - x) > rCur
                continue;
            end
            xCyl = x - cylCentre(1);
            z = sqrt(rCur^2 - xCyl^2);
            
            I(y,x) = z;
            
        end
    end
end

% I = imrotate(I,-30);
% slant = 0.1:0.1:30;
% Iadd = repmat(slant, [imSize,1]);
% I = I +Iadd;

%% create a mesh from this range image

[p,q] = meshgrid(1:imSize, 1:imSize);
ys = p(:);
xs = q(:);
ids = sub2ind(size(I), ys, xs);
zs = I(ids);

V = [xs,ys,zs];

lenF = 2*(imSize-1)^2;
lenV = size(V,1);
F = zeros(lenF,3);

cur = 1;

% imSize = 5;  % just to test triangulation
% lenV = 25;
% lenF = 2*(imSize-1)^2;
% F = zeros(lenF,3);

for i = 1:lenV
    
    if (i < lenV - imSize) && (mod(i, imSize) ~= 0)
        curFace = [i, i+1, i+imSize];
        F(cur,:) = curFace; 
        cur = cur + 1;
    end
    if (i > imSize) && (mod(i, imSize) ~= 1)
        curFace = [i, i-1, i - imSize];
        F(cur,:) = curFace; 
        cur = cur + 1;
    end
end


%% compute all parameters from the range image

radRC = 10;
minDepth = 1;
[r,c] = size(I);

halfImage(1) = floor(r/2)+1;
halfImage(2) = floor(c/2)+1;

%% visualize the results

% visualization parameters for the local frame of reference
vecLen = 25;
vectColors = [1,0,0; 0,1,0; 0,0,1]; 


% visualize the mesh
trisurf(F,V(:,1),V(:,2),V(:,3), 'FaceColor', [0.1, 0.1, 0.1], 'EdgeColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.4);
light('Position',[-1.0,-1.0,100.0],'Style','infinite');
axis equal;
lighting phong;
hold on

type = 'combinatorial';
options.normalize = 0;
W = compute_mesh_weight(V,F,type,options);



% compute geometry for this mesh
options.curvature_smoothing = 10;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(V,F,options);

inds = randperm(lenV);
subsetVis = 50;
inds = inds(1:subsetVis);

vectColors = eye(3);
vecLen = 30;

for i = 1:subsetVis 

    curPoint = V(inds(i), :);
    x = curPoint(1);
    y = curPoint(2);
    z = curPoint(3);
    
    % plot the eigenvectors in 3D
    
    VV = [Umin(:, inds(i))'; Umax(:, inds(i))'; Normal(:, inds(i))'];
    
    for ii = 1:3
        curVect = VV(ii,:);
        curColor = vectColors(ii, :);
        XX = [x, x + vecLen * curVect(1)];
        YY = [y, y + vecLen * curVect(2)];
        ZZ = [z, z + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor);
        hold on

    end
end
    
a = 2;
clear all



    


    
    
    
    
    
    
    
  


    

