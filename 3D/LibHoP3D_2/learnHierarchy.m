% This is a script to learn the entire 3D compositional hierarchy

% dataSetNumber 
% for Aim@Shape datasetnumber = 1;
% for Washington datasetnumber = 2

% partSelectionMethod
% for LibHob-motivated partSelectionMethod = 1 
% for optimization-based partSelectionMethod = 2


function [] = learnHierarchy()

dataSetNumber = 2;
partSelectionMethod = 1;
nClusters = 9;

% matlabpool open 4 % for parallel computing

if partSelectionMethod == 1
    abstractionLevel = 1;
elseif partSelectionMethod == 2
    
    % parameters for the different layers
    alphaParam = [0, 0, 0, 0];
    betaParam =  [0, 0, 1.5, 1.0];
    gammaParam = [0, 0, 0.33, 1.0];
end

% define folders configuration
root = '/home/vvk201/LibHoP3D/';

addpath([root,'utility']);
addpath([root,'settings']);
addpath([root,'Optimization']);
addpath([root,'Learning']);
addpath([root,'visualization']);
addpath([root,'statistics']);
addpath([root,'Recognition']);


% define the input data
if dataSetNumber == 1
    %depthPath = 'D:\3D\Input Data\Depth mapAllP\';
    depthPathDefault = [root, 'settings\list_depth.mat'];
    depthPath = '/home/vvk201/1TTT';
    
elseif dataSetNumber == 2
    depthPath_W = '/home/vvk201/Wash-rgbd-dataset_0003TT';
    depthPathDefault = '';
end


% output file names
%--------------------------------------------------------------------------
% files with co-occurrence statistics should have the following format
% statistics/statistics_layer_dataSetNumber_nClusters.mat

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);

if partSelectionMethod == 1
    aL = num2str(abstractionLevel);
end

% files for the raw statistics
statistics1Layer = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];
statistics3Layer = [root, 'statistics/statistics_3_', dsN, '_', nCl, '.mat'];
statistics4Layer = [root, 'statistics/statistics_4_', dsN, '_', nCl, '.mat'];
statistics5Layer = [root, 'statistics/statistics_5_', dsN, '_', nCl, '.mat'];
statistics6Layer = [root, 'statistics/statistics_6_', dsN, '_', nCl, '.mat'];

% files for the sieved statistics
statistics1LayerSieved = [root, 'statistics/statisticsSieved_1_', dsN, '_', nCl, '.mat'];
statistics3LayerSieved = [root, 'statistics/statisticsSieved_3_', dsN, '_', nCl, '.mat'];
statistics4LayerSieved = [root, 'statistics/statisticsSieved_4_', dsN, '_', nCl, '.mat'];
statistics5LayerSieved = [root, 'statistics/statisticsSieved_5_', dsN, '_', nCl, '.mat'];
statistics6LayerSieved = [root, 'statistics/statisticsSieved_6_', dsN, '_', nCl, '.mat'];

% files for the Aggregated statistics
statistics1LayerAggregated = [root, 'statistics/statisticsAggregated_1_', dsN, '_', nCl, '.mat'];
statistics3LayerAggregated = [root, 'statistics/statisticsAggregated_3_', dsN, '_', nCl, '.mat'];
statistics4LayerAggregated = [root, 'statistics/statisticsAggregated_4_', dsN, '_', nCl, '.mat'];
statistics5LayerAggregated = [root, 'statistics/statisticsAggregated_5_', dsN, '_', nCl, '.mat'];
statistics6LayerAggregated = [root, 'statistics/statisticsAggregated_6_', dsN, '_', nCl, '.mat'];

% files for part selection results
parts3Layer = [root, 'statistics/partsSelectionResults_3_', dsN, '_', nCl, '_a', aL, '.mat'];
parts4Layer = [root, 'statistics/partsSelectionResults_4_', dsN, '_', nCl, '_a', aL, '.mat'];
parts5Layer = [root, 'statistics/partsSelectionResults_5_', dsN, '_', nCl, '_a', aL, '.mat'];
parts6Layer = [root, 'statistics/partsSelectionResults_6_', dsN, '_', nCl, '_a', aL, '.mat'];

% files for overall coverage
coverageOverall3Layer = [root, 'statistics/coverageOverall_3_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverall4Layer = [root, 'statistics/coverageOverall_4_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverall5Layer = [root, 'statistics/coverageOverall_5_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverall6Layer = [root, 'statistics/coverageOverall_6_', dsN, '_', nCl, '_a', aL, '.mat'];

% files for the abstraction table
abstractionTable3Layer = [root, 'statistics/abstractionTable_3_', dsN, '_', nCl, '_a', aL, '.mat'];
abstractionTable4Layer = [root, 'statistics/abstractionTable_4_', dsN, '_', nCl, '_a', aL, '.mat'];
abstractionTable5Layer = [root, 'statistics/abstractionTable_5_', dsN, '_', nCl, '_a', aL, '.mat'];
abstractionTable6Layer = [root, 'statistics/abstractionTable_6_', dsN, '_', nCl, '_a', aL, '.mat'];


%---------define all parameters here----------------------------------------

fileListPrecomputed = false;
is_subset = false; % whether we shall use all files for learning

% define the subset length
if dataSetNumber == 1
    subset_len = 400; % how much shall we use for training
    subsetPercent = 1.0; % not used
elseif dataSetNumber == 2
    subset_len = 400;
    subsetPercent = 1.0; % what percent from each folder to use
end

% define parameters sigma and sigmaKernelSize

if dataSetNumber == 1
    ss = load('settings/sigma.mat');
    sigma = ss.sigma; % this is a sigma for gaussian derivative 
    sigmaKernelSize = ss.sigmaKernelSize; % size of the kernel for gaussian
elseif dataSetNumber == 2
    sigma = 2.0;
    sigmaKernelSize = 9;
end  

isErrosion = false;
discSize = 3;
dxKernel = load('dxKernel.mat');
dxKernel = dxKernel.dxKernel;
n2Clusters = nClusters^2;
combs = load('settings/combs10.mat');
combs = combs.combs; % combinations for the line discretization function
largestLine = 10; % for decomposition of the line discretization function

% ------------ parameters for line discretization---------------------
% min [error + wOverlap * overlap - wCoverage*coverage ]
if dataSetNumber == 1
    wCoverage = 0.25;
    wOverlap = 0.4;
elseif dataSetNumber == 2    
    wCoverage = 3.25;
    wOverlap = 0.52;
end

% --------Define what we learn here------------------------------------------
is_first_layer = false; % computes cluster centres, thresh and clusterSizes
is_third_layer = true;  % learns the first layer
is_4th_layer = false;
is_5th_layer = false; 
is_6th_layer = false; 

is_statistics3_collected = true;
is_statistics4_collected = false;
is_statistics5_collected = false;
is_statistics6_collected = false;

is_statistics3_sieved = true;
is_statistics4_sieved = false;
is_statistics5_sieved = false;
is_statistics6_sieved = false;

is_statistics3_aggregated = true;
is_statistics4_aggregated = true;
is_statistics5_aggregated = true;
is_statistics6_aggregated = true;


% here we define the first layer quantiles (bins of the first layer)
if nClusters == 9
    if dataSetNumber == 1
        quantilesFirst = [0.045 0.14 0.27 0.38];
    elseif dataSetNumber == 2
        quantilesFirst = [0.045 0.12 0.24 0.32];
    end
elseif nClusters == 7
    if dataSetNumber == 1
        quantilesFirst = [0.05 0.17 0.33];
    elseif dataSetNumber == 2
        quantilesFirst = [0.045 0.19 0.30];
    end
end

%--------------------------------------------------------------------------
% creating a filelist here

if dataSetNumber == 1
    [list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);
    list_mask = [];
elseif dataSetNumber == 2
    [list_depth, list_mask, ~, lenF] = extractFileListWashington(fileListPrecomputed, depthPath_W, depthPathDefault, is_subset, subsetPercent);
end

% % % %  downsampling (if necessary)
% if dataSetNumber == 1
%     is_downsampling = false;
%     dowsample_rate = 1;
% elseif dataSetNumber == 2    
%     is_downsampling = true;
%     dowsample_rate = 2.5;
% end
% 
% upsampleImages(list_depth, list_mask, lenF, is_downsampling, dowsample_rate); % to be done only once!

is_downsampling = false;
dowsample_rate = 1;

% fieldSize = [15, 5, 71];  % x, y and z directions
% coverageOverall = ProjectAllStatistics3(statistics3LayerSieved, statistics3LayerAggregated, fieldSize, list_depth, lenF);
% save(coverageOverall3Layer,  'coverageOverall');

%--------------------------------------------------------------------------

if is_first_layer  % here we learn the first layer parameters
    disp('Calibration of the first layer parameters...');

    [cluster1Centres, cluster1Lengths, thresh] = learnFirstLayer(list_depth, list_mask, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, ...
                                            nClusters, is_downsampling, dowsample_rate, quantilesFirst, dataSetNumber);

    save(statistics1Layer,  'cluster1Centres', 'cluster1Lengths', 'thresh', 'nClusters', 'dataSetNumber');
end

%--------------------------------------------------------------------------
if is_third_layer  % here we lear the third layer of the hierarchy
    
    disp('Learning of the third layer ...');
    % we have to perform the following procedures
    
    if ~is_first_layer
        % read the first layer
        load(statistics1Layer);
    end
    
    %load('displacements_Layers3_4.mat');
    load('settings/displacements_Layer3.mat');
    lenDisp = size(displacements, 1);   
    fieldSize = [17, 5, 71];  % x, y and z directions
    depthStep = thresh/3;
    quant = 0.06;
    
    
    if ~is_statistics3_collected
        [outputStatistics, outputCoords, curTS] = CollectStats_3_Layer(list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discSize, nClusters, ...
                                                        cluster1Centres, cluster1Lengths, thresh, combs, largestLine, displacements, lenDisp, ... 
                                                        wCoverage, wOverlap, fieldSize, depthStep, is_downsampling, dowsample_rate, dataSetNumber);

        save(statistics3Layer, 'outputStatistics', 'outputCoords');
    end
    
    if ~is_statistics3_sieved 
    
        load(statistics3Layer);
        curTS = size(outputStatistics, 1);
        

        statistics = outputStatistics;
        clear('outputStatistics');
        curTS = size(statistics, 1);
        lenDisp = 2;

        % for spead up reasons we sort statistics by the first column (central element)
        [~, order] = sort(statistics(:,1));
        statistics = statistics(order,:);
        outputCoords = outputCoords(order,:);
        clear('order');

        % compute the most frequent pairs
        thresh3Pair = 0.02 * lenF;
        [tablePairs, numPairs] = Learn3Pairs(statistics, curTS, n2Clusters, thresh3Pair, lenDisp);

        % filter out rows with the least frequent pairs. those pairs with less than thresh3Pair occurrances will be filtered out
        [ind, statistics, curTS] = Sieve3Pairs(statistics, curTS, tablePairs, lenDisp);
        outputCoords = outputCoords(ind, :);

        % now we have to measure depth and eliminate rows with depth discontinuities
        [cluster3Depths] = compute3Depths(statistics, n2Clusters, quant, lenDisp);

        % now we have to filter out the depths with discontinuity
        [ind, statistics, curTS] = Sieve3PairsD(statistics, curTS, cluster3Depths, lenDisp);
        outputCoords = outputCoords(ind,:);
        
        save(statistics3LayerSieved, 'statistics', 'cluster3Depths', 'outputCoords');
        
    else
        
        % read the sieved statistics
%         load('statistics/statistics3Layer_4000_sieved.mat');

        load(statistics3LayerSieved);  %  variable statistics should be read
        curTS = size(statistics, 1);
    end
    
    % aggregate statistics and prepare it to be fed to the optimization
    % function
    
    if ~is_statistics3_aggregated % aggregate statistics
        
        sieve_thresh = 3;
        [X, frequencies, curTS, triples] = Aggregate_3Layer(statistics, nClusters, n2Clusters, curTS, sieve_thresh);
        
        % sieve 'statistics' and 'outputCoords' once again
        [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieve_thresh, X);
        outputCoords = outputCoords(ind, :);  
        
        save(statistics3LayerSieved, 'statistics', 'cluster3Depths', 'outputCoords');
        save(statistics3LayerAggregated, 'X' ,'frequencies', 'curTS', 'triples');
    end
    
    if partSelectionMethod == 1
       
        checkImages(list_depth, lenF); % makes Images in the folder 3 channels ones
        
        [partsOut, coverageOut, lenOut] = partSelectionNew(nClusters, n2Clusters, statistics3LayerSieved, statistics3LayerAggregated, ...
                                 dataSetNumber, fieldSize, list_depth, lenF, abstractionLevel, abstractionTable3Layer);
        save(parts3Layer, 'partsOut', 'coverageOut', 'lenOut', 'abstractionLevel');

    elseif partSelectionMethod == 2
        
        load(statistics3LayerSieved);    
        load(statistics3LayerAggregated);
    
        [triples3Out] = Optimization_layer3(X, frequencies, triples, nClusters, alphaParam(3), betaParam(3), gammaParam(3), subset_len); % this is the main procedure
        save('statistics/layer3.mat', 'triples3Out');
                
        n3Clusters = size(triples3Out, 1);

%       [triple3OutDepth] = store3Layer(triples3Out, cluster3Depths, n3Clusters, nClusters);
%       save('statistics/layer3.mat', 'triple3OutDepth', 'triples3Out');

    end
    
    % store the visualized vocabulary to this folder
    
    %str_folder = 'D:\3D\Demonstrations\layer3Elements\';
    %[is_ok] = layer_3_demonstrator(triple3OutDepth, nClusters, n3Clusters, str_folder, fieldSize, depthStep, cluster1Centres);
    
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
%         load('statistics/Layer3_final.mat');    % read triples3Out
%         load('statistics/triple3OutDepth.mat'); % triple3OutDepth

        load('statistics/layer3_Wash1.mat');  % triples3Out
        n3Clusters = size(triples3Out, 1);
    end
    
    load('displacements_Layers3_4.mat');
    lenDisp = size(displacements, 1);   
    fieldSize = [13, 13, 71];
    depthStep = thresh/3;
    quant = 0.07;
    maxDist = 2;
    
  [outputStatistics, curTS] = CollectStats_3_4_Layers(list_depth, lenF, sigma, sigmaKernelSize, dxKernel, nClusters, cluster1Centres, cluster1Lengths, thresh, ...
   outRoot, combs, largestLine, displacements, wCoverage, wOverlap, fieldSize, depthStep);
    
    is_sieved = true;

    
    %----------------------------------------------------------------------
    % HERE WE LEARN THE FOURTH LAYER
    
    % first we sieve the outputStatistics
    if ~is_sieved
        load('statistics/statistics3Layer_Wash_Overall.mat'); % outputStatistics
        [outputStatistics, curTS] = Sieve4LayerStatistics(outputStatistics, lenF, n2Clusters, quant); 
    else
       load('statistics/statistics4Layer_sieved_Wash'); % outputStatistics
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
    thresh4Pair = 0.002 * lenF;
    [tablePairs, numPairs] = Learn3Pairs(statistics, curTS, n3Clusters, thresh4Pair, lenDisp);
    
    % filter out rows with the least frequent pairs. those pairs with less than thresh3Pair occurrances will be filtered out
    [~, statistics, curTS] = Sieve3Pairs(statistics, curTS, tablePairs, lenDisp);
    
    % now we have to measure depth and eliminate rows with depth discontinuities
    [cluster3Depths] = compute3Depths(statistics, n3Clusters, quant, lenDisp);
    
    % now we have to filter out the depths with discontinuity
    [~, statistics, curTS] = Sieve3PairsD(statistics, curTS, cluster3Depths, lenDisp);
    
    sieve_thresh = 0;
    [X, frequencies, curTS, triples4] = Aggregate_4Layer(statistics, n3Clusters, triples3Out, curTS, sieve_thresh);   
    
    [triples4Out] = Optimization_layer4(X, frequencies, triples3Out, triples4, nClusters, n3Clusters, maxDist, subset_len, alphaParam(4), betaParam(4), gammaParam(4));
    n4Clusters = size(triples4Out, 1);
    
    save('statistics/layer4Wash1.mat', 'triples4Out');
    a = 2;
      
    
end   



% matlabpool close 

