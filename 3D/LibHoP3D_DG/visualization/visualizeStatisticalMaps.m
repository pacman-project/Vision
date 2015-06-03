function [] = visualizeStatisticalMaps( layerID, tripleOutDepth, displ, ...
                                            nClusters, nNClusters, str_folder, fieldSize, depthStep, cluster1Centres, isFIG, stats5D)

    
    fieldCenter = ceil(fieldSize / 2);
    f = figure('Position', [100, 100, 700, 900]);
            
    for jjj = 1:nNClusters% for each 6th layer element
        

        [positions, elements] = partMeanReconstruction(layerID, jjj, fieldCenter, tripleOutDepth, displ{3}, displ{5}, displ{7}, nClusters);
         
        for jj = 1:nNClusters
            
            % put a central part in the middle of field
            subplot(3,1,1:2)
            surfaceVisualizerT(fieldSize, positions, elements, nClusters, cluster1Centres, depthStep, [0,0,1], false, 0.3);
            hold on
            
            
            curStats = stats5D{jjj}{jj};
            
            maxV = max(max(max(curStats)));
            if maxV == 0 
                hold off;
                clf(f);
                continue;
            end
            
            if isempty(curStats)
                hold off;
                clf(f);
                continue;
            end

            % XYZS are variables used for 3D scatter plot (C - color, S - size)
            % first initialize central point
            X = [1, fieldCenter(1), fieldSize(1)]; Y = [1, fieldCenter(2), fieldSize(2)]; Z = [1, fieldCenter(3), fieldSize(3)]; 
            S = [1,5,1];
            C = rand(3);

            colormap cool;
            cmap = colormap; % cmap contains 64 colors

            for i = 1:fieldSize(1)
                for j = 1:fieldSize(2)
                    for k = 1:fieldSize(3)
                        curValue = curStats(i,j,k);
                        if curValue > 0 
                            
                            if curValue == maxV
                                indX = i;
                                indY = j;
                                indZ = k;  % position with the highest frequency
                            end
                            X = [X, i];
                            Y = [Y, j];
                            Z = [Z, k];
                            S = [S, ceil(60 * (curValue)/(maxV)+1)]; % size of the current marker
                            ind = ceil(64 * (curValue)/(maxV));
                            C = [C; cmap(ind,:)];
                        end
                    end
                end
            end
            Z = Z*depthStep;
            scatter3(X,Y,Z,S,C, 'filled');
                        
            subplot(3,1,3)
%             % now visualize the other part
            [positionsAdd, elementsAdd] = partMeanReconstruction(layerID, jj, fieldCenter, tripleOutDepth, displ{3}, displ{5}, displ{7}, nClusters);
%             % shift to the right position
%             positionsAdd(1) = indX;
%             positionsAdd(2) = indY;
%             positionsAdd(3) = indZ;
%             
            surfaceVisualizerT(fieldSize, positionsAdd, elementsAdd, nClusters, cluster1Centres, depthStep, [1,0,0], false, 0.2);
            
            
            A = exist(str_folder, 'dir');
            if A == 0
                mkdir(str_folder);
            end

            if ~ isFIG
                str = ['layer_', num2str(layerID), '_', num2str(jjj), '_', num2str(jj), '.png'];
                str1 = [str_folder, str];
                saveas(f, str1, 'png');
            else
                str = ['layer_', num2str(layerID), '_', num2str(jjj), '_', num2str(jj), '.fig'];
                str1 = [str_folder, str];
                saveas(f, str1, 'fig');
            end
            
            hold off;
            clf(f);
            
            % may also apply "rotate" function for better observation
            
        end
        


    end
    
    is_ok = true;
end


