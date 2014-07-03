% this is to display layer 4 elements one after another

function [is_ok] = layer_5_demonstrator(triple5OutDepth, triple4OutDepth, triple3OutDepth, displ3, displ5, nClusters, n5Clusters, ...
                                        str_folder, fieldSize, depthStep, cluster1Centres)

 
    offsetY4 = [-displ3, 0, displ3];
    offsetX5 = [-displ5, 0, displ5];

    % pre-compute a table
    table = zeros(nClusters, nClusters);

    for i = 1:nClusters
        for j = 1:nClusters
            table(i,j) = compute2elementIndex(i, j, nClusters);
        end
    end
    
    indsEls4 = [1,5,6];

    f = figure;
            
    for j = 1:n5Clusters  % for each 5th layer element
        
        cur5triple = triple5OutDepth(j, :);
        els5 = cur5triple(indsEls4);
        depths5Adder = [cur5triple(4), 0 cur5triple(9)];  % top centre bottom
        
        
        for k = 1:3 % three triples of the 4th layer
        
            cur4triple = triple4OutDepth(els5(k), :);
            els4 = cur4triple(indsEls4);
            depths4Adder = [cur4triple(4), 0 cur4triple(9)];  % top centre bottom
        
            for i = 1:3                            % three triples of the third layer
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

                positionLeft   = [fieldCenter(1) - displ3 + offsetX5(k), fieldCenter(2) + offsetY4(i), fieldCenter(3) + depths(3) + depths4Adder(i) + depths5Adder(k)];
                positionCentre = [fieldCenter(1)          + offsetX5(k), fieldCenter(2) + offsetY4(i), fieldCenter(3)             + depths4Adder(i) + depths5Adder(k)];
                positionRight =  [fieldCenter(1) + displ3 + offsetX5(k), fieldCenter(2) + offsetY4(i), fieldCenter(3) + depths(6) + depths4Adder(i) + depths5Adder(k)];

                positions = [positionLeft; positionCentre; positionRight]; % left, middle, right
                surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);
                hold on

            end
            
        end
        
        
        
        A = exist(str_folder, 'dir');

        if A == 0
            mkdir(str_folder);
        end

        str = ['layer5_', num2str(j), '.png'];
        str1 = [str_folder, str];
        saveas(f, str1, 'png');

        hold off

    %         str = ['layer3_', num2str(i), '.fig'];
    %         str1 = [str_folder, str];
    %         savefig(str1);
    end
    
    is_ok = true;
