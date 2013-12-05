% this applies ISODATA to the Third layer!!!
stats = load('../statistics/Layer3Aggregated.mat');
stats = stats.triples;

% rewrite to the format required for ISODATA
X = zeros(10000000,6);
frequencies = zeros(1,10000000);
% embedding of X into 6 dimensional space

thresh = 5;

% pre-computing table

table = zeros(2,81);
for i = 1:81
    [clusterX, clusterY] = compute2derivatives(i, 9);
    table(1,i) = clusterX;
    table(2,i) = clusterY;
end

ind = 0;
for i = 1:81
    for j = 1:81
        for k = 1:81
            freq = stats(i,j,k);
            if freq > thresh
                ind = ind + 1;
                X(ind,1) = table(1, i);
                X(ind,2) = table(2, i); % left
                X(ind,3) = table(1, j);
                X(ind,4) = table(2, j); % center
                X(ind,5) = table(1, k);
                X(ind,6) = table(2, k); % right
                frequencies(ind) = freq; 
            end
        end
    end
end

disp('data converted');
frequencies = frequencies(1:ind);
X = X(1:ind, :);
rangeMin = 1;
rangeMax = 9;
kInit = 200;
nMin =  3000; % minimal number of points to form the cluster
iterations = 10;
sigmaMax = 1.0;
LMin = 4;
pMax = 3;
want_disp = false;

[Z] = Isodata_Vladislav(X, frequencies, rangeMin, rangeMax ,kInit, nMin, iterations, sigmaMax, LMin, pMax, want_disp, 3);

% Final = ISODATA(X', kInit, kInit, nMin, sigmaMax, LMin, pMax, iterations);

