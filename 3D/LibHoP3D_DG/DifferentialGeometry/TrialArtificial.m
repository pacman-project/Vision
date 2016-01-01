% this is to try different differential geometry algorthms for the
% artificial data
function TrialArtificial()

elementType = 2;
elementRadius = 15;
is2D = false;

imSize = 300;
I = zeros(imSize, imSize);

methodId = 2; 

isIteration = 1;
numIter = 2;
%% here we create a primitive
primType = 1; % cylinder
zeroThresh = 10^(-3);
slant = pi/6;
elevation = 1;

if primType == 1 % cylinder
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
    
elseif primType == 3 % Planar
    
    r = 70;
    h = 200;
    z0 = 10;
    
    cylCentre = [150, 150, 0];
    cylinderStartX = cylCentre(1) - r;
    cylinderEndX = cylCentre(1) + r;
    cylinderStartY = cylCentre(2) - h/2;
    cylinderEndY = cylCentre(2) + h/2;
    
    for y = cylinderStartY:cylinderEndY
        for x = cylinderStartX:cylinderEndX
            
            
            z = z0 + x * 0.5;
            
            I(y,x) = z;
            
        end
    end
    
elseif primType == 4  % primitive is uploaded from the file
    
    I = imread('D:\Input Data\VladislavSTD\Vladislav_STD\DataSet\depth_0.50\99_5_2_1_3_2.png');
    sc = load('D:\LibHoP3D_DG\settings\calibration3.mat');
    zScale = sc.zScale;
    I = I(:,:,1);
    I = double(I);
    I = I*zScale;
    
elseif primType == 5 % trancated cone
    
    rt = 10;
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
% imtool(Iadd, [1,200]);

% I(I<22) = 20; 


%% here we compute all the transformations



%% compute derivatives (if it is required)
imtool(I, [min(min(I)), max(max(I))]);

gf = fspecial('gaussian', 5 ,1.0);
Ig = imfilter(I, gf);
dx = 0.5*[-1 0 1];
dy = dx';

Ix = imfilter(Ig, dx);
Iy = imfilter(Ig, dy);


%% compute all parameters from the range image

radRC = 10;
minDepth = 1;
[r,c] = size(I);

halfImage(1) = floor(r/2)+1;
halfImage(2) = floor(c/2)+1;

%% visualize the results

% visualization parameters for the local frame of reference
vecLen = 40;
vectColors = [1,0,0; 0,1,0; 0,0,1]; 

if primType == 4
    I = I-400;
    I(I<0) = 0;
%     I = flipud(I);
end


surf(I, 'FaceColor',[0.3, 0.3 0.3], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'FaceLighting', 'phong');
camlight left
axis equal;
hold on


 for x = 150:150        % 80:10:220% 90:10:210
    for y = 150:150     % 100:40:180

        curDepth = I(y,x);
        [indsXOut, indsYOut, depths] = computeNeighbors(I, x, y, elementRadius, is2D, minDepth);
        
    %   test for a neihbourhood found 
        scatter3(indsXOut, indsYOut, depths);


    % convert these coordinates to the local frame of
    %  reference (centered at the point [j,i,curDepth])
        xs = indsXOut - x;
        ys = indsYOut - y;
        zs = depths - I(y,x);
        
        xs = xs(2:end);
        ys = ys(2:end);
        zs = zs(2:end);
        
%         imtool(I, [min(min(I)), max(max(I))]);

%         ids1 = find(xs < 0);
%         ids2 = find( xs>=0);
%         perm = randperm(length(ids1));
%         perm = perm(1:length(ids1)/5);
%         xs = [xs(ids2), xs(perm)];
%         ys = [ys(ids2), ys(perm)];
%         zs = [zs(ids2), zs(perm)];
%         figure; scatter3(xs, ys, zs); axis equal;
        

        if methodId == 1  % paraboloid fitting

            [H, K, QF1, QF2, S, V, D, is_ok] = paraboloidFitting(xs, ys, zs);                        
            plotFrame(V, vecLen, vectColors, x,y,curDepth);

        elseif methodId == 2 % method of Gabriel Taubin (modified)

            [V, D] = computeDarbouxFrame([0,0,1]', xs, ys, zs); 
            plotFrame(V, vecLen, vectColors, x,y,curDepth);

        end

            % plot the eigenvectors in 3D

        a = 2;
    
    end 
end
 
end

function plotFrame(V, vecLen, vectColors, x,y,curDepth)
    for ii = 1:3
        curVect = V(:, ii);
        curColor = vectColors(ii, :);
        XX = [x, x + vecLen * curVect(1)];
        YY = [y, y + vecLen * curVect(2)];
        ZZ = [curDepth, curDepth + vecLen * curVect(3)];

        plot3(XX, YY, ZZ, 'Color', curColor);
        hold on

    end
end






% this is an attemp to make iterative approach

%             if isIteration
%                 for i = 1:numIter 
%                     % convert points to the computed frame of reference and
%                     % repeat computations
% 
%                     Xtemp = V(:,2);
%                     Ytemp = V(:,3);  % local x and y axis
%                     NormalCentral = V(:,1);
% 
%                     T = eye(4);
%                     R = eye(4,4); R(1:3, 1) = Xtemp'; R(1:3, 2) = Ytemp'; R(1:3, 3) = NormalCentral';
%                     M = inv(R)*T; 
% 
%                     vect = M * [xs;  ys;  zs;  ones(1, length(xs))];
% 
%                     % try again in the new frame of reference
%                     xs = vect(1,:);
%                     ys = vect(2,:);
%                     zs = vect(3,:);
%                     [H, K, QF1, QF2, S, V, D] = paraboloidFitting(xs, ys, zs);
% 
%                     % return V to the original (global) frame of reference
%                     V = [V; ones(1,3)];
%                     V = R * V;
%                     V = V(1:3, 1:3);
% 
%                     % plot the eigenvectors in 3D
%                     plotFrame(V, vecLen, vectColors, x,y,curDepth)
% 
%                     a = 2;
%                 end
%             end

    




    


    
    
    
    
    
    
    
  


    