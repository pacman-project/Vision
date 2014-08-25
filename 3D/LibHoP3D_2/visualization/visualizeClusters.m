clear all;


% stats = load('statistics5D_3layer_19_10000.mat');
% stats = stats.statistics;

cluster3stats = load('statistics/pairs3layer.mat');
cluster3stats = cluster3stats.cluster3stats;
set(0,'defaultfigureposition',[100, 200, 560 420]);


gridX = 19;
gridY = 19;
gridZ = 71;

lenC =size(cluster3stats,1);

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
% [clusterCenters, clusterSize] = defineCluster1Centers(nClusters, thresh);

% % to see what are class probabilities
% sumC = zeros(nClusters, nClusters);
% 
% for i = 1:nClusters
%     for j = 1:nClusters
%         clusterN = (i - 1) * nClusters + j;
%         sumC(i,j) = sum(sum(sum(sum(sum(stats(clusterN,:,:,:,:)))))); 
%     end
% end

for x1 = 5:5
    for y1 = 5:5
        for x2 = 5:5
            for y2 = 4:6

                clusterX1 = x1;
                clusterY1 = y1;
                clusterX2 = x2;
                clusterY2 = y2;

                el1 = (clusterX1 - 1) * nClusters + clusterY1;
                el2 = (clusterX2 - 1) * nClusters + clusterY2;

                % we have to find rows with first element el1 and second el2 (up to 4 raws) 
                miniStats = [];
                for i = 1:lenC
                    if cluster3stats(i,1) == el1 && cluster3stats(i,2) == el2
                        miniStats = [miniStats; cluster3stats(i,:)]; 
                    end
                end

                lenC1 = size(miniStats,1);


                [x,y] = meshgrid(8:1:12, 8:1:12);

                dx1 = 0.8 * clusterCenters(clusterX1)/depthStep;
                dy1 = 0.7 * clusterCenters(clusterY1)/depthStep;

                z = dx1 * (x) + dy1 * y;
                z = z - z(3,3);
                z = z + ceil(gridZ/2);
                figure

                mesh(x,y,z, 'FaceColor','blue','EdgeColor','none', 'FaceAlpha', 1);
                axis([1, gridX, 1, gridY, 1, gridZ]);
                hold on

                len = size(cluster3stats,1); % amount of clusters
                % take the first one just to see what is going on
                for i = 1:lenC1
                    mu = miniStats(i, 3:5);
                    S = miniStats(i, 6:14);
                    Sigma = reshape(S, 3, 3);
                    drawCluster(mu, Sigma);
                    hold on
                end

                hold on

                mesh(x,y,z, 'FaceColor','blue','EdgeColor','none', 'FaceAlpha', 1);
                axis([1, gridX, 1, gridY, 1, gridZ]);

                hold off;
                
            end
        end
    end
end
