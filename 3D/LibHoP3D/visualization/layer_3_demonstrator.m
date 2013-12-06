% this is to display layer 3 elements one after another

function [is_ok] = layer_3_demonstrator(triple3OutDepth, nClusters, n3Clusters, str_folder, fieldSize, depthStep, cluster1Centres)

    if (nargin < 1)
        Z = load('statistics/triple3OutDepth.mat');
        triple3OutDepth = Z.triple3OutDepth;
        nClusters = 9;
        n3Clusters = size(triple3OutDepth, 1);
        str_folder = 'D:\3D\Demonstrations\layer3Elements\';
        fieldSize = [15, 15, 71];
        thresh = 112;
        depthStep = thresh / 5;
    end

    % pre-compute a table
    table = zeros(nClusters, nClusters);

    for i = 1:nClusters
        for j = 1:nClusters
            table(i,j) = compute2elementIndex(i, j, nClusters);
        end
    end

    f = figure;

    for i = 1:n3Clusters                             % for each third layer element
        curEl = triple3OutDepth(i,:);
        indsEl = [1,2,6,7,8,9];
        indsDepths = [3,4,5,10,11,12];
        
        curXY = curEl(indsEl);
        depths = curEl(indsDepths);
        
        elements = zeros(1,3);
        for j = 1:3
            elements(j) = table(curXY(2*j-1), curXY(2*j));
        end
        % define positions
        fieldCenter = ceil(fieldSize / 2);
        positionLeft = [fieldCenter(1) - 4, fieldCenter(2), fieldCenter(3) + depths(3)];
        positionRight =[fieldCenter(1) + 4, fieldCenter(2), fieldCenter(3) + depths(6)];
        
        positions = [positionLeft; fieldCenter; positionRight]; % left, middle, right
        surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);

        str = ['layer3_', num2str(i), '.png'];
        str1 = [str_folder, str];
        saveas(f, str1);

    end
    
    is_ok = true;
        

