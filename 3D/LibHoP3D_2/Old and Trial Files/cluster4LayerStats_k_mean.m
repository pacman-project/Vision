clear all;

scaleFactor = 10;

% this is to aggregate the 4th layer statistics
% by means of k-mean clustering
nClusters = 9;
statistics = load('statistics/Statistics5D_4layer_19_2300');
stats = statistics.tripleStatistics;
len = size(stats,1); % number of samples
len = round (len / scaleFactor);

stats = stats(1:len, :);

features = zeros(len, 18);

matrixXY = zeros(nClusters, nClusters);
matrix2 = zeros(nClusters*nClusters, 2);

for i = 1:nClusters
    for j = 1:nClusters
        ind = compute2elementIndex(i, j, nClusters);
        matrixXY(i,j) = ind;
        matrix2(ind,1) = i;
        matrix2(ind,2) = j;
    end
end  

% here we rewrite our feature matrix
for i = 1:len
    for j = 1:9
        curEl = stats(i,j);
        clusterX = matrix2(curEl, 1);
        clusterY = matrix2(curEl, 2);
        features(i, 2*j - 1: 2*j) = [clusterX, clusterY]; 
    end
end

% now we perform clustering of this statistics to n4Clusters
n4Clusters = 500;
[IDX,C,sumD] = kmeans(features, n4Clusters, 'emptyaction', 'drop', 'Display', 'iter', 'maxIter', 50);


avgDist = sum(sumD) / len;






