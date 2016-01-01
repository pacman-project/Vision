% this applies ISODATA to the Third layer (used for edge and corner features)

function [clusterCentres, sampleLables] = ISODATA_3_Layer(tableParts, nClusters)

    n2Clusters = nClusters^2;
    lenTP= size(tableParts,1);

    % data used for clustering
    X = zeros(lenTP,4);
    frequencies = zeros(lenTP, 1);

    SieveThresh = 5;

    % pre-computing table

    table = zeros(2,n2Clusters);
    for i = 1:n2Clusters
        [clusterX, clusterY] = compute2derivatives(i, nClusters);
        table(1,i) = clusterX;
        table(2,i) = clusterY;
    end

    for i = 1:lenTP
        freq = tableParts(i,3);
        nI = tableParts(i,4);  % number of images where this part appear

        if freq > SieveThresh
            X(i,1) = table(1, tableParts(i,1));
            X(i,2) = table(2, tableParts(i,1)); % left
            X(i,3) = table(1, tableParts(i,2));
            X(i,4) = table(2, tableParts(i,2)); % right
            frequencies(i) = freq; 
        end
    end

    rangeMin = 1;
    rangeMax = 9;
    kInit = 500;
    nMin =  100; % minimal number of points to form the cluster
    iterations = 100;
    sigmaMax = 1.0;
    LMin = 2;
    pMax = 10;
    want_disp = false;
    is_discrete = true;

    [clusterCentres, sampleLables] = Isodata_Vladislav(X, frequencies, rangeMin, rangeMax ,kInit, nMin, iterations, sigmaMax, LMin, pMax, want_disp, is_discrete);

end

