% clear all;

stats = load('../Statistics/Statistics5D_4layer_19_2300_aggregated.mat');
stats = stats.tripleNew;

X = stats(:,1:9);
[r,c] = size(X);
XX = zeros(r,18);
% embedding of X into 18 dimensional space

% pre-computing table
table = zeros(3,81);
table(1,:) = 1:81;
for i = 1:81
    [clusterX, clusterY] = compute2derivatives(i, 9);
    table(2,i) = clusterX;
    table(3,i) = clusterY;
end

for i = 1:r
    curLine = X(i,:);
    for j = 1:9
        el = X(i,j);
        XX(i, 2*j - 1) = table(2, el);
        XX(i, 2*j)     = table(3, el);
    end
end

disp('data converted');
clear('X', 'table');

frequencies = stats(:,10);
rangeMin = 1;
rangeMax = 9;
kInit = 200;
nMin =  1000; % minimal number of points to form the cluster
iterations = 20;
sigmaMax = 1.2;
LMin = 7;
pMax = 30;
want_disp = false;

[Z] = Isodata_Vladislav(XX, frequencies, rangeMin, rangeMax ,kInit, nMin, iterations, sigmaMax, LMin, pMax, want_disp);

[r,c] = size(Z);

