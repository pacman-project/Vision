% this function is purposed to learn the parameters of the first layer
% applicable for the Washington dataset only

function [cluster1Centres, cluster1Bounds, thresh] = learnFirstLayer(list_depth, list_mask, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, ...
                      nClusters, quantilesFirst, dataSetNumber, is_guided, r_guided, eps, ...
                      is_mask_extended, maxExtThresh1, maxExtThresh2)
    
    if ~((nClusters == 9)||(nClusters == 7))
        disp('Error. Please define quantiles for this number of clusters');
    end
    
    allDXs = zeros(10^7,1); % here we collect all the values of dx
    
    if dataSetNumber == 1 || dataSetNumber == 3
        [allDXs, lenX] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize,...
                              is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2);
    elseif dataSetNumber == 2
        [allDXs, lenX] = computeFirstLayerDistribution_W(list_depth, list_mask, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize,...
                              is_guided, r_guided, eps, is_mask_extended, maxExtThresh1, maxExtThresh2);
    end
    
    allDXs = allDXs(1:lenX);
    
    if nClusters == 9 || nClusters == 7
        [bounds] = quantile(allDXs, quantilesFirst);
    end
    
    %% Here we process statistics defining cluster centres and lengths
    
    cluster1Bounds = [bounds, fliplr(abs(bounds))]; % to make it symmetric (if for whatever reason it is not)
    cluster1Centres = zeros(1, nClusters);
    for i = 1:nClusters
        cluster1Centres(i) = (cluster1Bounds(i) + cluster1Bounds(i+1)) / 2;
    end
    
    cluster1Centres(1) = cluster1Bounds(2) - 0.01;
    cluster1Centres(nClusters) = cluster1Bounds(nClusters) + 0.01;
    thresh = -cluster1Bounds(2);
    
end



%     cluster1Centres = [cluster1Centres, 0, fliplr(abs(cluster1Centres))]; % to make it symmetric (if for whatever reason it is not)
% %    cluster1Centres = round(cluster1Centres);
% 
%     thresh = cluster1Centres(end);
%     
%     % solve a small system to compute cluster lengths here
%     
%     X = [];
%     a = zeros(nClusters,1);
%     for i = 1:nClusters - 1
%         curX = zeros(1,nClusters);
%         curX(i:i+1) = ones(1,2)/2;
%         X = [X; curX];
%         a(i) = cluster1Centres(i+1) - cluster1Centres(i);
%     end
%     
%     centInd = floor(nClusters/2);
%     centralInterval = cluster1Centres(centInd + 2) - cluster1Centres(centInd);
%     
%     % now the last condition
%     a(end) = centralInterval*0.4; % length of the central element    
%     curX = zeros(1,nClusters);
%     curX(ceil(nClusters/2)) = 1;
%     X = [X; curX]; 
%     cluster1Lengths = linsolve(X,a);