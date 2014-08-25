% this function is to visualize the raw statistical maps

clear all;


stats = load('statistics/statistics5D_3layer_19_100.mat');
stats = stats.statistics;
set(0,'defaultfigureposition',[-700, 200, 560 420]);

gridX = 19;
gridY = 19;
gridZ = 71;

centerX = ceil(gridX/2);
centerY = ceil(gridY/2);
centerZ = ceil(gridZ/2);

initialDepth = 26;

% it's a 81 * 81 * gridX * gridY * gridZ  variable

% first normalize the statistics
% sumAll = sum(sum(sum(sum(sum(stats)))));
% stats = stats/sumAll;
% sumAll = sum(sum(sum(sum(sum(stats)))));
% stats(stats<0.000001) = 0;


nClusters = 9;
thresh = 95;
n2Clusters = nClusters^2;
depthStep = thresh/3;
middleCluster = ceil(nClusters/2);
 

% cluster sizes in percent (if we have nClusters = 9)
if nClusters == 9
    clusterSizes = load('settings/cluster1Sizes.mat'); 
    clusterSizes = clusterSizes.clusterSizes;
    cluster1SizesPercent = clusterSizes* 0.01;
end

[clusterCenters, cluster1Length] = defineCluster1Centers(nClusters, cluster1SizesPercent, thresh);


for x1 = 5:5
    for y1 = 5:5
        for x2 = 5:5
            for y2 = 5:6

                clusterX1 = x1;
                clusterY1 = y1;
                clusterX2 = x2;
                clusterY2 = y2;
                

                el1 = (clusterX1 - 1) * nClusters + clusterY1;
                el2 = (clusterX2 - 1) * nClusters + clusterY2;
                
                curStats = stats(el1, el2, :, :, :);
                curStats =  squeeze(curStats);
                maxV = max(max(max(curStats)));
                if maxV == 0
                    continue;
                end
                
                figure
                
                % XYZS are variables used for 3D scatter plot (C - color, S - size)
                % first initialize central point
                X = [1, centerX, gridX]; Y = [1, centerY, gridY]; Z = [1, centerZ, gridZ]; S = [1,5,1];
                C = rand(3);

                colormap cool;
                cmap = colormap; % cmap contains 64 colors
                
                for i = 1:gridX
                    for j = 1:gridY
                        for k = 1:gridZ
                            curValue = curStats(i,j,k);
                            if curValue > 0 
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
                
                [x,y] = meshgrid(centerX-2:1:centerX+2, centerY-2:1:centerX+2);
                dx1 = 0.8 * clusterCenters(clusterX1)/depthStep;
                dy1 = 0.8 * clusterCenters(clusterY1)/depthStep;
                z = dx1 * (x) + dy1 * y;
                z = z - z(3,3);
                z = z + ceil(gridZ/2);

                mesh(x,y,z,'FaceColor','blue','EdgeColor','none', 'FaceAlpha', 0.5);
                axis([1, gridX, 1, gridY, 1, gridZ]);

                hold on
                scatter3(X,Y,Z,S,C, 'filled');
                % colorbar('location','southoutside');


                hold off;
            end
        end
    end
end 



%                 Cmax = max(C);
%                 Cmin = min(C);
%                 shiftC = Cmax - Cmin;
%                 colorStep = 63 / shiftC;
%                 S = floor((C - Cmin) * colorStep) + 1; % S - color of
%                 every marker
