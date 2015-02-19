% this is to display layer 4 elements one after another

function [is_ok] = layer_6_demonstrator(triple6OutDepth, triple5OutDepth, triple4OutDepth, triple3OutDepth, displ3, displ5, nClusters, n6Clusters, ...
                                        str_folder, fieldSize, depthStep, cluster1Centres)

    fieldCenter = ceil(fieldSize / 2);
    f = figure;
            
    for j = 1:n6Clusters  % for each 6th layer element
        
        [positions, elements] = partMeanReconstruction(6, j, fieldCenter, [], [], triple6OutDepth, triple5OutDepth, triple4OutDepth, ...
                                                    triple3OutDepth, displ3, displ5, 0, nClusters);
        surfaceVisualizerT(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);
        A = exist(str_folder, 'dir');

        if A == 0
            mkdir(str_folder);
        end

        str = ['layer6_', num2str(j), '.png'];
        str1 = [str_folder, str];
        saveas(f, str1, 'png');

%         str = ['layer3_', num2str(i), '.fig'];
%         str1 = [str_folder, str];
%         savefig(str1);
    end
    
    is_ok = true;
