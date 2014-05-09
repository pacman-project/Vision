% this function is purposed to learn the parameters of the first layer
% applicable for the Washington dataset only

function [cluster1Centres, cluster1Lengths, thresh] = learnFirstLayer(list_depth, list_mask, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, ...
                                            nClusters, is_downsampling, dowsample_rate, quantilesFirst, dataSetNumber)
    
    if ~((nClusters == 9)||(nClusters == 7))
        disp('Error. Please define quantiles for this number of clusters');
    end
    
    allDXs = zeros(10^7,1); % here we collect all the values of dx
    
    if dataSetNumber == 1
        [allDXs, lenX] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize,...
                                                                                            is_downsampling, dowsample_rate);
    elseif dataSetNumber == 2
        [allDXs, lenX] = computeFirstLayerDistribution_W(list_depth, list_mask, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize,...
                                                                                            is_downsampling, dowsample_rate);
    end
    
    allDXs = allDXs(1:lenX);
    
    if nClusters == 9 || nClusters == 7
        [cluster1Centres] = quantile(allDXs, quantilesFirst);
    end
    
    cluster1Centres = [cluster1Centres, 0, fliplr(abs(cluster1Centres))]; % to make it symmetric (if for whatever reason it is not)
%    cluster1Centres = round(cluster1Centres);

    thresh = cluster1Centres(end);
    
    % solve a small system to compute cluster lengths here
    
    X = [];
    a = zeros(nClusters,1);
    for i = 1:nClusters - 1
        curX = zeros(1,nClusters);
        curX(i:i+1) = ones(1,2)/2;
        X = [X; curX];
        a(i) = cluster1Centres(i+1) - cluster1Centres(i);
    end
    
    centInd = floor(nClusters/2);
    centralInterval = cluster1Centres(centInd + 2) - cluster1Centres(centInd);
    
    % now the last condition
    a(end) = centralInterval*0.3; % length of the central element    
    curX = zeros(1,nClusters);
    curX(ceil(nClusters/2)) = 1;
    X = [X; curX]; 
    cluster1Lengths = linsolve(X,a);

end