% this is to display layer 3 elements one after another

function [is_ok] = layer_3_demonstrator(triple3OutDepth, displ, nClusters, n3Clusters, str_folder, fieldSize, depthStep, cluster1Centres)


    f = figure;
    fieldCenter = ceil(fieldSize / 2);

    for i = 1: n3Clusters % for each third layer element
        
        [positions, elements] = partMeanReconstruction(3, i, fieldCenter, [], [], [], [], [], ...
                                                            triple3OutDepth, displ, 0, 0, nClusters);
                                                        
        surfaceVisualizerT(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);

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

        a = 2;
    end
    
    is_ok = true;
        

