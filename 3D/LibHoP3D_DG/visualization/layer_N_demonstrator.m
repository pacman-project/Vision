% this function demonstrates vocabularies of all layers

function [is_ok] = layer_N_demonstrator(layerID, tripleOutDepth, offsetsConventional, ...
                                            nClusters, nNClusters, str_folder, fieldSize, depthStep, cluster1Centres, isFIG)

    fieldCenter = ceil(fieldSize / 2);
    f = figure;
            
    for j = 1:nNClusters  % for each 6th layer element
        
        [positions, elements] = partMeanReconstruction(layerID, j, fieldCenter, tripleOutDepth, offsetsConventional, depthStep, nClusters);
                                                       
        
        surfaceVisualizerT(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep);
        A = exist(str_folder, 'dir');

        if A == 0
            mkdir(str_folder);
        end

        if ~ isFIG
            str = ['layer', num2str(layerID), '_', num2str(j), '.png'];
            str1 = [str_folder, str];
            saveas(f, str1, 'png');
        else
            str = ['layer', num2str(i), '.fig'];
            str1 = [str_folder, str];
            saveas(f, str1, 'fig');
        end
    end
    
    is_ok = true;
end