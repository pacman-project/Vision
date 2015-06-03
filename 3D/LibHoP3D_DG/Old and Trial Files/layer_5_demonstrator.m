% this is to display layer 4 elements one after another

function [is_ok] = layer_5_demonstrator(triple5OutDepth, triple4OutDepth, triple3OutDepth, displ3, displ5, nClusters, n5Clusters, ...
                                        str_folder, fieldSize, depthStep, cluster1Centres)

    fieldCenter = ceil(fieldSize / 2);
    f = figure;
            
    for j = 1:n5Clusters  % for each 5th layer element
        
        [positions, elements] = partMeanReconstruction(5, j, fieldCenter, [], [], [], triple5OutDepth, triple4OutDepth, ...
                                                    triple3OutDepth, displ3, displ5, 0, nClusters);
        surfaceVisualizerT(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);
        A = exist(str_folder, 'dir');

        if A == 0
            mkdir(str_folder);
        end

        str = ['layer5_', num2str(j), '.png'];
        str1 = [str_folder, str];
        saveas(f, str1, 'png');

%         str = ['layer3_', num2str(i), '.fig'];
%         str1 = [str_folder, str];
%         savefig(str1);
    end
    
    is_ok = true;
