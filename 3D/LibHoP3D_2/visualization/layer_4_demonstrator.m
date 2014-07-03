% this is to display layer 4 elements one after another

function [is_ok] = layer_4_demonstrator(triple4OutDepth, triple3OutDepth, displ, nClusters, n4Clusters, str_folder, fieldSize, depthStep, cluster1Centres)

    if (nargin < 1)
%         Z = load('statistics/triple4OutDepth.mat');
%         triple3OutDepth = Z.triple4OutDepth;
%         nClusters = 9;
%         displ = 6;
%         n3Clusters = size(triple3OutDepth, 1);
%         str_folder = 'D:/3D/Demonstrations/layer3Elements/';
%         fieldSize = [17, 17, 71];
%         thresh = 98;
%         depthStep = thresh / 3;
    end
    
    offsetY = [-displ, 0, displ];

    % pre-compute a table
    table = zeros(nClusters, nClusters);

    for i = 1:nClusters
        for j = 1:nClusters
            table(i,j) = compute2elementIndex(i, j, nClusters);
        end
    end
    
    indsEls4 = [1,5,6];

    f = figure;
            
    for j = 1:n4Clusters  % for each 4th layer element
        
        
        cur4triple = triple4OutDepth(j, :);
        els4 = cur4triple(indsEls4);
        depths4Adder = [cur4triple(4), 0 cur4triple(9)];  % top centre bottom
        
        for i = 1:3                            % three triples of the third layet
            curEl = triple3OutDepth(els4(i),:);
            indsEl = [1,2,6,7,8,9];
            indsDepths = [3,4,5,10,11,12];

            curXY = curEl(indsEl);
            depths = curEl(indsDepths);

            elements = zeros(1,3);
            for jj = 1:3
                elements(jj) = table(curXY(2*jj-1), curXY(2*jj));  % part index in range [0, n2Clusters]
            end
            % define positions
            fieldCenter = ceil(fieldSize / 2);
            
            positionLeft   = [fieldCenter(1) - displ, fieldCenter(2) + offsetY(i), fieldCenter(3) + depths(3) + depths4Adder(i)];
            positionCentre = [fieldCenter(1)        , fieldCenter(2) + offsetY(i), fieldCenter(3)             + depths4Adder(i)];
            positionRight =  [fieldCenter(1) + displ, fieldCenter(2) + offsetY(i), fieldCenter(3) + depths(6) + depths4Adder(i)];

            positions = [positionLeft; positionCentre; positionRight]; % left, middle, right
            surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);
            hold on

        end
        
        A = exist(str_folder, 'dir');

        if A == 0
            mkdir(str_folder);
        end

        str = ['layer4_', num2str(j), '.png'];
        str1 = [str_folder, str];
        saveas(f, str1, 'png');

        hold off

    %         str = ['layer3_', num2str(i), '.fig'];
    %         str1 = [str_folder, str];
    %         savefig(str1);
    end
    
    is_ok = true;
