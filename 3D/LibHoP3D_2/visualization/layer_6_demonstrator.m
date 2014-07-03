% this is to display layer 4 elements one after another

function [is_ok] = layer_6_demonstrator(triple6OutDepth, triple5OutDepth, triple4OutDepth, triple3OutDepth, displ3, displ5, nClusters, n6Clusters, ...
                                        str_folder, fieldSize, depthStep, cluster1Centres)

 
    offsetY4 = [-displ3, 0, displ3];
    offsetX5 = [-displ5, 0, displ5];
    offsetY6 = [-displ5, 0, displ5];

    % pre-compute a table
    table = zeros(nClusters, nClusters);

    for i = 1:nClusters
        for j = 1:nClusters
            table(i,j) = compute2elementIndex(i, j, nClusters);
        end
    end
    
    indsEls4 = [1,5,6];

    f = figure;
            
    for j = 1:n6Clusters  % for each 6th layer element
        
        cur6triple = triple6OutDepth(j, :);
        els6 = cur6triple(indsEls4);
        depths6Adder = [cur6triple(4), 0 cur6triple(9)];  % top centre bottom
        
        
        for kk = 1:3 % three triples of the layer 5
            
            cur5triple = triple5OutDepth(els6(kk), :);
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

                    positionLeft   = [fieldCenter(1) - displ3 + offsetX5(k), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3) + depths(3) + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk)];
                    positionCentre = [fieldCenter(1)          + offsetX5(k), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3)             + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk)];
                    positionRight =  [fieldCenter(1) + displ3 + offsetX5(k), fieldCenter(2) + offsetY4(i) + offsetY6(kk), fieldCenter(3) + depths(6) + depths4Adder(i) + depths5Adder(k) + depths6Adder(kk)];

                    positions = [positionLeft; positionCentre; positionRight]; % left, middle, right
                    surfaceVisualizer(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);
                    hold on

                end

            end
        end
        
        
        
        A = exist(str_folder, 'dir');

        if A == 0
            mkdir(str_folder);
        end

        str = ['layer6_', num2str(j), '.png'];
        str1 = [str_folder, str];
        saveas(f, str1, 'png');

        hold off

    %         str = ['layer3_', num2str(i), '.fig'];
    %         str1 = [str_folder, str];
    %         savefig(str1);
    end
    
    is_ok = true;
