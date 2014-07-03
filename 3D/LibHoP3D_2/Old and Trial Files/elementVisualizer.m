triple3stats = load('triplesNew.mat'); % load('tirpleStatistics5D_3layer_19_5000.mat');
triple3stats = triple3stats.tripleNew;

cluster3stats = load('cluster3stats.mat');
cluster3stats = cluster3stats.cluster3stats;

firstColumn = cluster3stats(:,1);
secondColumn = cluster3stats(:,2);
lastColumn = cluster3stats(:,16);


% positions: [x1,y1,z1; x2,y2,z2; ...] may be real numbers far from margin
% elements (35,44,51,...)
% nCluster - integer
% thresh - real number (refers to thresh form the 1st layer)
% fieldSize is a vector [sizeX, SizeY, sizeZ]
% depthStep - real number (presumably 95/5)
% function [out] = surfaceVisualizer(fieldSize, positions, elements, nClusters, thresh, depthStep)

    nClusters = 9;
    thresh = 95;
    depthStep = thresh/3;
    fieldSize = [19,19,71];
    positions = zeros(5,3);
    positions(1,:) = [10,10,36];

for i = 40:50
    line = triple3stats(i,:);
    elements = line(1:5);
    el = line(1);
    
    score = line(6);
    if score > 200
    end
    % we have 5 elements totally
    % now we define positions
    for j = 2:5
        curEl = line(j);
        indices = find(firstColumn == el & secondColumn == curEl& lastColumn == j);
        indices = indices(1);
        bigLine = cluster3stats(indices, :);
        mu = bigLine(3:5);
        positions(j,:) = mu;
    end
    
    figure
    [out] = surfaceVisualizer(fieldSize, positions, elements, nClusters, thresh, depthStep);
end