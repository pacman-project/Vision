% This is a script to learn the entire 3D compositional hierarchy

% define folders configuration
root = 'D:\3D\LibHoP3D\';

addpath([root,'utility']);
addpath([root,'settings']);
addpath([root,'Optimization']);
addpath([root,'Learning']);
addpath([root,'visualization']);

% matlabpool open 4 % for parallel computing


% define the input data
depthPath = 'D:\3D\Input Data\Depth mapAllP\';
depthPathDefault = [root, 'settings\list_depth.mat'];
fileListPrecomputed = true; 
is_subset = true; % whether we shall use all files for learning
subset_len = 4000; % how much shall we use


% define and download the parameters
ss = load('settings/sigma.mat');
sigma = ss.sigma; % this is a sigma for gaussian derivative 
sigmaKernelSize = ss.sigmaKernelSize; % size of the kernel for gaussian
dxKernel = load('dxKernel.mat');
dxKernel = dxKernel.dxKernel;
nClusters = 9;
n2Clusters = nClusters^2;
combs = load('settings/combs10.mat');
combs = combs.combs; % combinations for the line discretization function
largestLine = 10; % for decomposition of the line discretization function

% ------------ parameters for line discretization---------------------
% min [error + wOverlap * overlap - wCoverage*coverage ]
wCoverage = 0.25;
wOverlap = 0.4;

% --------Define what we learn here------------------------------------------
is_first_layer = false; % computes cluster centres, thresh and clusterSizes
is_third_layer = false;  % learns the first layer
is_4th_layer = true;
is_5th_layer = false; 
is_6th_layer = false; 

%--------------------------------------------------------------------------
% forming a filelist here

[list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);


%--------------------------------------------------------------------------

if is_first_layer  % here we learn the first layer parameters
    disp('Calibration of the first layer parameters...');
    [cluster1Centres, cluster1Lengths, thresh] = learnFirstLayer(list_depth, sigma, sigmaKernelSize, dxKernel, nClusters);
    %[allDXs, lenX] = computeFirstLayerDistribution(list_depth, sigma, sigmaKernelSize, dxKernel);   
end

%--------------------------------------------------------------------------
if is_third_layer  % here we lear the third layer of the hierarchy
    
    disp('Learning of the third layer ...');
    % we have to perform the following procedures
    
    if ~is_first_layer
        % read the first layer
        struct = load([root, 'settings/firstLayer']);
        cluster1Centres = struct.cluster1Centres;
        cluster1Lengths = struct.cluster1Lengths;
        thresh = struct.thresh;
    end
    
    load('displacements_Layers3_4.mat');
    lenDisp = size(displacements, 1);   
    fieldSize = [13, 13, 71];
    depthStep = thresh/4;
    quant = 0.10;
    
    
%     [outputStatistics, curTS] = CollectStats_3_4_Layers(list_depth, lenF, sigma, sigmaKernelSize, dxKernel, nClusters, cluster1Centres, cluster1Lengths, thresh, ...
%                  outRoot, combs, largestLine, displacements, wCoverage, wOverlap, fieldSize, depthStep);
    
    load('statistics/statistics3Layer_4000.mat');
    curTS = size(outputStatistics, 1);
   
    
    %----------------------------------------------------------------------
    % HERE WE LEARN THE THIRD LAYER
    
    % displacements = [0,0; 0,-4; 4,0; 0,4; -4,0; 4,-4; 4,4; -4,4; -4,-4];
    
    %      9 5 8      
    %      2 1 4      
    %      6 3 7
    
    a = [9,5,8; 2,1,4; 6,3,7];
    
    statistics = [];
    % prepare statistics for the third layer (use all three rows together)
    for i = 1:3
        cur = a(i,:);   % ex.  9,5,8
        cols = [(cur(2)-1)*2, (cur(1)-1)*2, (cur(1)-1)*2+1, (cur(3)-1)*2, (cur(3)-1)*2+1];
        cols(cols == 0) = 1; % if i == 2 (middle element)
        curStatistics = outputStatistics(:, cols);
        if i ~= 2       % we have to recompute relative depths
            absDepthsCol = (cur(2)-1)*2 + 1;
            absDepths = outputStatistics(:, absDepthsCol);
            curStatistics(:,3) = curStatistics(:,3) - absDepths;
            curStatistics(:,5) = curStatistics(:,5) - absDepths;
        end
        statistics = [statistics; curStatistics];
    end
    clear('outputStatistics');
    curTS = size(statistics, 1);
    lenDisp = 2;
    
    % for spead up reasons we sort statistics by the first column (central element)
    [~, order] = sort(statistics(:,1));
    statistics = statistics(order,:);
    clear('order');

    % compute the most frequent pairs
    thresh3Pair = 0.02 * lenF;
    [tablePairs, numPairs] = Learn3Pairs(statistics, curTS, n2Clusters, thresh3Pair, lenDisp);
    
    % filter out rows with the least frequent pairs. those pairs with less than thresh3Pair occurrances will be filtered out
    [~, statistics, curTS] = Sieve3Pairs(statistics, curTS, tablePairs, lenDisp);
    
    % now we have to measure depth and eliminate rows with depth discontinuities
    [cluster3Depths] = compute3Depths(statistics, n2Clusters, quant, lenDisp);
    
    % now we have to filter out the depths with discontinuity
    [~, statistics, curTS] = Sieve3PairsD(statistics, curTS, cluster3Depths, lenDisp);
    
    % now we come to learning of triples!!! (in X-direction)
%     load('statistics/statistics3Layer_3300_filtered.mat');
%     curTS = size(statistics, 1);
    
    % aggregate statistics and prepare it to be fed to the optimization
    % function
    
    sieve_thresh = 1;
    [X, frequencies, curTS, triples] = Aggregate_3Layer(statistics, nClusters, n2Clusters, curTS, sieve_thresh);
    [triples3Out] = Optimization_layer3(X, frequencies, triples, nClusters); % this is the main procedure
    
    n3Clusters = size(triples3Out, 1);
    
    [triple3OutDepth] = store3Layer(triples3Out, cluster3Depths, n3Clusters, nClusters);
    
    str_folder = 'D:\3D\Demonstrations\layer3Elements\';
    [is_ok] = layer_3_demonstrator(triple3OutDepth, nClusters, n3Clusters, str_folder, fieldSize, depthStep, cluster1Centres);
    
    % here we can also save the results
end



if is_4th_layer
    
    disp('Learning of the 4th layer ...');
    
    if ~is_first_layer
        % read the first layer
        struct = load([root, 'settings/firstLayer']);
        cluster1Centres = struct.cluster1Centres;
        cluster1Lengths = struct.cluster1Lengths;
        thresh = struct.thresh;
    end
    
    if ~is_third_layer
        load('statistics/Layer3_final.mat'); %read triples3Out
        load('statistics/triple3OutDepth.mat'); % triple3OutDepth
        n3Clusters = size(triples3Out, 1);
    end
    
    load('displacements_Layers3_4.mat');
    lenDisp = size(displacements, 1);   
    fieldSize = [13, 13, 71];
    depthStep = thresh/3;
    quant = 0.07;
    numDisp = 6;
    
%   [outputStatistics, curTS] = CollectStats_3_4_Layers(list_depth, lenF, sigma, sigmaKernelSize, dxKernel, nClusters, cluster1Centres, cluster1Lengths, thresh, ...
%    outRoot, combs, largestLine, displacements, wCoverage, wOverlap, fieldSize, depthStep);
    
    is_sieved = true;

    
    %----------------------------------------------------------------------
    % HERE WE LEARN THE FOURTH LAYER
    
    % first we sieve the outputStatistics
    if ~is_sieved
        load('statistics/statistics3Layer_4000.mat'); % outputStatistics
        [outputStatistics, curTS] = Sieve4LayerStatistics(outputStatistics, lenF, n2Clusters, quant); 
    else
       load('statistics/statistics3Layer_4000_sieved.mat');
       curTS = size(outputStatistics, 1);
    end
    
    outputStatistics = int16(outputStatistics);
    %  Replace horisontal triples with the third layer elements
    %      9 5 8      
    %      2 1 4      
    %      6 3 7  --->   [central, top, top_depth, bottom, bottom_depth]
    
    [statistics, curTS] = replaceTripleWith3Elements(outputStatistics, curTS, triples3Out, n2Clusters, nClusters, true);
    lenDisp = 2;
    
    % for spead up reasons we sort statistics by the first column (central element)
    [~, order] = sort(statistics(:,1));
    statistics = statistics(order,:);
    clear('order');
    
    % compute the most frequent pairs
    thresh4Pair = 0.003 * lenF;
    [tablePairs, numPairs] = Learn3Pairs(statistics, curTS, n3Clusters, thresh4Pair, lenDisp);
    
    % filter out rows with the least frequent pairs. those pairs with less than thresh3Pair occurrances will be filtered out
    [~, statistics, curTS] = Sieve3Pairs(statistics, curTS, tablePairs, lenDisp);
    
    % now we have to measure depth and eliminate rows with depth discontinuities
    [cluster3Depths] = compute3Depths(statistics, n3Clusters, quant, lenDisp);
    
    % now we have to filter out the depths with discontinuity
    [~, statistics, curTS] = Sieve3PairsD(statistics, curTS, cluster3Depths, lenDisp);
    
    sieve_thresh = 0;
    [X, frequencies, curTS, triples4] = Aggregate_4Layer(statistics, n3Clusters, triples3Out, curTS, sieve_thresh);
    [triples4Out] = Optimization_layer4(X, frequencies, triples3Out, triples4, nClusters, n3Clusters, numDisp, subset_len);
    
    a = 2;
      
    
end


    
    
    



% matlabpool close 

