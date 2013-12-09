% this function is purposed to learn the parameters of the first layer

function [cluster1Centres, cluster1Lengths, thresh] = learnFirstLayer(list_depth, sigma, sigmaKernelSize, dxKernel, nClusters)
    
    if nClusters ~= 9
        disp('Error. Please define quantiles for this number of clusters');
    end
    
    allDXs = zeros(10^7,1); % here we collect all the values of dx
    %allDYs = zeros(10^7,1);
    
    %[allDXs, allDYs, lenX, LenY] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel); 
    
    [allDXs, lenX] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel);
    allDXs = allDXs(1:lenX);
    %allDYs = allDYs(1:lenY);
    
    if nClusters == 9
        [cluster1Centres] = quantile(allDXs,[0.045 0.14 0.28 0.40]);
    end
    
    cluster1Centres = [cluster1Centres, 0, fliplr(abs(cluster1Centres))]; % to make it symmetric (if for whatever reason it is not)
    [cluster1Centres] = round(cluster1Centres);
    thresh = cluster1Centres(end);
    
    % solve a small system to compute cluster lengths here
    
    X = [];
    a = zeros(nClusters,1);
    for i = 1:nClusters - 1
        curX = zeros(1,9);
        curX(i:i+1) = ones(1,2)/2;
        X = [X; curX];
        a(i) = cluster1Centres(i+1) - cluster1Centres(i);
    end
    % now the last condition
    a(end) = thresh*0.07; % length of the central element    
    curX = zeros(1,9);
    curX(ceil(nClusters/2)) = 1;
    X = [X; curX]; 
    cluster1Lengths = linsolve(X,a);

end