% this is to display layer 3 elements one after another

function [is_ok] = layer_3_demonstrator(triple3OutDepth, displ, nClusters, n3Clusters, str_folder, fieldSize, depthStep, cluster1Centres)

    if (nargin < 1)
        Z = load('statistics/triple3OutDepth.mat');
        triple3OutDepth = Z.triple3OutDepth;
        nClusters = 9;
        displ = 6;
        n3Clusters = size(triple3OutDepth, 1);
        str_folder = 'D:/3D/Demonstrations/layer3Elements/';
        fieldSize = [17, 17, 71];
        thresh = 98;
        depthStep = thresh / 3;
    end

    % pre-compute a table
    table = zeros(nClusters, nClusters);

    for i = 1:nClusters
        for j = 1:nClusters
            table(i,j) = compute2elementIndex(i, j, nClusters);
        end
    end

    f = figure;

    for i = 1: n3Clusters                             % for each third layer element
        curEl = triple3OutDepth(i,:);
        indsEl = [1,2,6,7,8,9];
        indsDepths = [3,4,5,10,11,12];
        
        curXY = curEl(indsEl);
        depths = curEl(indsDepths);
        
        elements = zeros(1,3);
        for j = 1:3
            elements(j) = table(curXY(2*j-1), curXY(2*j));  % part index in range [0, n2Clusters]
        end
        % define positions
        fieldCenter = ceil(fieldSize / 2);
        positionLeft = [fieldCenter(1) - displ, fieldCenter(2), fieldCenter(3) + depths(3)];
        positionRight =[fieldCenter(1) + displ, fieldCenter(2), fieldCenter(3) + depths(6)];
        
        positions = [positionLeft; fieldCenter; positionRight]; % left, middle, right
        surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);

        A = exist(str_folder, 'dir');
    
        if A == 0
            mkdir(str_folder);
        end
        
        str = ['layer3_', num2str(i), '.png'];
        str1 = [str_folder, str];
        saveas(f, str1, 'png');
        
%         str = ['layer3_', num2str(i), '.fig'];
%         str1 = [str_folder, str];
%         savefig(str1);


    end
    
    is_ok = true;
        

