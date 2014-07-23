% this is to display layer 4 elements one after another

function [is_ok] = layer_4_demonstrator(triple4OutDepth, triple3OutDepth, displ, nClusters, n4Clusters, str_folder, fieldSize, depthStep, cluster1Centres)

    f = figure;
    fieldCenter = ceil(fieldSize / 2);
            
    for j = 1:n4Clusters  % for each 4th layer element
        
        
        [positions, elements] = partMeanReconstruction(4, j, fieldCenter, [], [], [], [], triple4OutDepth, ...
                                                            triple3OutDepth, displ, 0, 0, nClusters);
        surfaceVisualizerT(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);

        
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
