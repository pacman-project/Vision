function [ is_ok ] = VisualizeCornerParts(cluster1Centres)

    faceColor = [0,0,1];
    transparency = 0.9;
    edgeTransparency = transparency;
    isFlip = false;
    fieldSize = [17, 17, 17];
    centre = ceil(fieldSize/2);
    str_folder = 'D:\Visualized vocabulary\Vladislav_STD\3_layer_loc\';
    isFIG = false;
    layerID = 3;
    partIDStart = 50;

    f = figure;
    for i = 5:5%1:length(cluster1Centres)
        
        curSlope = pi *(90 + cluster1Centres(i)) / 180;
        
        
        len = 5;
        point1 = [6 + len * cos(curSlope), 1, 1 + len * sin(curSlope)];
        point2 = [6 + len * cos(curSlope), 6, 1 + len * sin(curSlope)];
        
        
        V = [1,1,1;  1,6,1;  6,6,1; 6,1,1; point1; point2];
        
        V = rotateVerts(degtorad(0), degtorad(0), degtorad(90), V);
        
        
        for j = 1:size(V, 1)
            V(j,:) = V(j,:) + centre - 3;
        end
        
        F = [1,2,3; 3,4,1; 5,6,4; 3,4,6];
        
        
        trisurf(F, V(:,1),V(:,2),V(:,3),'FaceColor',faceColor, 'EdgeColor', faceColor, 'FaceAlpha', transparency, 'EdgeAlpha', edgeTransparency);
        xlabel('x')
        ylabel('y')
        axis equal; 
        axis([0, fieldSize(1), 0, fieldSize(2), 0, fieldSize(3), 50, 60]); % ceil(1 + scale * rangeZ), ceil(rangeZ - scale * rangeZ)
        if isFlip
            set(gca,'zdir','reverse')
        end

        A = exist(str_folder, 'dir');

        if A == 0
            mkdir(str_folder);
        end

        if ~ isFIG
            str = ['layer', num2str(layerID), '_', num2str(i + partIDStart), '.png'];
            str1 = [str_folder, str];
            saveas(f, str1, 'png');
        else
            str = ['layer', num2str(layerID), '_', num2str(i + partIDStart), '.fig'];
            str1 = [str_folder, str];
            saveas(f, str1, 'fig');
        end
    end
    
    is_ok = true;
        


end

function radAlpha = degtorad(alpha)
    radAlpha = pi*(alpha)/180;
end

function V1 = rotateVerts(angleX, angleY, angleZ, V1)

    rotZ = [cos(angleZ), -sin(angleZ), 0,0; sin(angleZ), cos(angleZ), 0,0; 0,0,1,0; 0,0,0,1];
    rotX = [1,0,0,0; 0 cos(angleX), -sin(angleX), 0; 0, sin(angleX), cos(angleX),0; 0,0,0,1];
    rotY = [cos(angleY), 0, sin(angleY), 0; 0,1,0,0; -sin(angleY), 0,cos(angleY), 0; 0,0,0,1];
    
    numPoints = size(V1, 1);
    R = rotX * rotY * rotZ;
    
    for jj = 1:numPoints
        Vcur = [V1(jj,:), 1];
        Vcur = R*Vcur';
        V1(jj, :) = Vcur(1:3)';
    end
  
end

