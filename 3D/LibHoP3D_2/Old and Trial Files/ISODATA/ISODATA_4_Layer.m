% clear all;

Z = load('../Statistics/layer3.mat');
Z = Z.Z;  % these are 3rd layer elements

[rZ, cZ] = size(Z);

% this applies ISODATA to the fourth layer!!!
stats = load('../statistics/Layer4Aggregated.mat');
stats = stats.triples;

% rewrite to the format required for ISODATA
X = zeros(10000000,3);
frequencies = zeros(1,10000000);
thresh = 5;

ind = 0;
for i = 1:rZ
    for j = 1:rZ
        for k = 1:rZ
            freq = stats(i,j,k);
            if freq > thresh
                ind = ind + 1;
                X(ind,:) = [i,j,k];
                frequencies(ind) = freq; 
            end
        end
    end
end

disp('data converted');
frequencies = frequencies(1:ind);
X = X(1:ind, :);
rangeMin = 1;
rangeMax = rZ;
kInit = 400;
nMin =  3000; % minimal number of points to form the cluster
iterations = 40;
sigmaMax = 1.0;
LMin = 9;
pMax = 4;
want_disp = false;

[Z] = Isodata_Vladislav(X, frequencies, rangeMin, rangeMax ,kInit, nMin, iterations, sigmaMax, LMin, pMax, want_disp, 4);

