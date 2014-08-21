% This is a script to learn the entire 3D compositional hierarchy

% dataSetNumber 
% for Aim@Shape dataSetNumber = 1;
% for Washington dataSetNumber = 2;
% for Vladislav's_STD dataSetNumber = 3;
% grasping dataset dataSetNumber = 4;


% partSelectionMethod
% for LibHob-motivated partSelectionMethod = 1 
% for optimization-based partSelectionMethod = 2



function [] = learnHierarchy()

dataSetNumber = 2;
dataSetNames{1} = 'Aim@Shape';
dataSetNames{2} = 'Washington';
dataSetNames{3} = 'Vladislav_STD';

partSelectionMethod = 1;
nClusters = 7;

n2Clusters = nClusters^2;

% here we initialize the parallel computing (including GPU acceleration)
% g = gpuDevice(1);
% reset(g);
% matlabpool open 10    % for parallel computing

if partSelectionMethod == 1
    
    [ meargeThresh3, meargeThresh4,meargeThresh5, meargeThresh6, meargeThresh7, meargeThresh8 ] = defineMeargeThreshes(1);
    
elseif partSelectionMethod == 2
    
    % parameters for the different layers
    alphaParam = [0, 0, 0, 0];
    betaParam =  [0, 0, 1.5, 1.0];
    gammaParam = [0, 0, 0.33, 1.0];
end

% define folders configuration
commonRoot = 'D:/'; 
root = [commonRoot, 'LibHoP3D/'];
addPaths(root);

depthPathDefault = '';


%% define the input data
if dataSetNumber == 1
    depthPath = [commonRoot, 'Input Data/AimShape/4T_600'];   %'D:\3D\Input Data\Images for categorization\1T';     
elseif dataSetNumber == 2
    depthPath = [commonRoot, 'Input Data/Washington/Washington3Categories_008'];  %'/home/vvk201/Wash-rgbd-dataset'; 
elseif dataSetNumber == 3
    depthPath = [commonRoot, 'commonRoot/Input Data/VladislavSTD/Vladislav_STD_600/depth'];
end


%% define output file names

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);

if partSelectionMethod == 1
    aL = '3';
end

% files for the raw statistics
statistics1Layer = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];
statistics3Layer = [root, 'statistics/statistics_3_', dsN, '_', nCl, '.mat'];
statistics4Layer = [root, 'statistics/statistics_4_', dsN, '_', nCl, '.mat'];
statistics5Layer = [root, 'statistics/statistics_5_', dsN, '_', nCl, '.mat'];
statistics6Layer = [root, 'statistics/statistics_6_', dsN, '_', nCl, '.mat'];
statistics7Layer = [root, 'statistics/statistics_7_', dsN, '_', nCl, '.mat'];
statistics8Layer = [root, 'statistics/statistics_8_', dsN, '_', nCl, '.mat'];

% files for the sieved statistics
statistics1LayerSieved = [root, 'statistics/statisticsSieved_1_', dsN, '_', nCl, '.mat'];
statistics3LayerSieved = [root, 'statistics/statisticsSieved_3_', dsN, '_', nCl, '.mat'];
statistics4LayerSieved = [root, 'statistics/statisticsSieved_4_', dsN, '_', nCl, '.mat'];
statistics5LayerSieved = [root, 'statistics/statisticsSieved_5_', dsN, '_', nCl, '.mat'];
statistics6LayerSieved = [root, 'statistics/statisticsSieved_6_', dsN, '_', nCl, '.mat'];
statistics7LayerSieved = [root, 'statistics/statisticsSieved_7_', dsN, '_', nCl, '.mat'];
statistics8LayerSieved = [root, 'statistics/statisticsSieved_8_', dsN, '_', nCl, '.mat'];

% files for the Aggregated statistics
statistics1LayerAggregated = [root, 'statistics/statisticsAggregated_1_', dsN, '_', nCl, '.mat'];
statistics3LayerAggregated = [root, 'statistics/statisticsAggregated_3_', dsN, '_', nCl, '.mat'];
statistics4LayerAggregated = [root, 'statistics/statisticsAggregated_4_', dsN, '_', nCl, '.mat'];
statistics5LayerAggregated = [root, 'statistics/statisticsAggregated_5_', dsN, '_', nCl, '.mat'];
statistics6LayerAggregated = [root, 'statistics/statisticsAggregated_6_', dsN, '_', nCl, '.mat'];
statistics7LayerAggregated = [root, 'statistics/statisticsAggregated_7_', dsN, '_', nCl, '.mat'];
statistics8LayerAggregated = [root, 'statistics/statisticsAggregated_8_', dsN, '_', nCl, '.mat'];

% files for part selection results
parts3Layer = [root, 'statistics/partsSelectionResults_3_', dsN, '_', nCl, '_a', aL, '.mat'];
parts4Layer = [root, 'statistics/partsSelectionResults_4_', dsN, '_', nCl, '_a', aL, '.mat'];
parts5Layer = [root, 'statistics/partsSelectionResults_5_', dsN, '_', nCl, '_a', aL, '.mat'];
parts6Layer = [root, 'statistics/partsSelectionResults_6_', dsN, '_', nCl, '_a', aL, '.mat'];
parts7Layer = [root, 'statistics/partsSelectionResults_7_', dsN, '_', nCl, '_a', aL, '.mat'];
parts8Layer = [root, 'statistics/partsSelectionResults_8_', dsN, '_', nCl, '_a', aL, '.mat'];

% files for overall coverage
coverageOverall3Layer = [root, 'statistics/coverageOverall_3_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverall4Layer = [root, 'statistics/coverageOverall_4_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverall5Layer = [root, 'statistics/coverageOverall_5_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverall6Layer = [root, 'statistics/coverageOverall_6_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverall7Layer = [root, 'statistics/coverageOverall_7_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverall8Layer = [root, 'statistics/coverageOverall_8_', dsN, '_', nCl, '_a', aL, '.mat'];

% files for layer visualization
fileForVisualization3Layer = [root, 'statistics/fileForVisualization_3_', dsN, '_', nCl, '.mat'];
fileForVisualization4Layer = [root, 'statistics/fileForVisualization_4_', dsN, '_', nCl, '.mat'];
fileForVisualization5Layer = [root, 'statistics/fileForVisualization_5_', dsN, '_', nCl, '.mat'];
fileForVisualization6Layer = [root, 'statistics/fileForVisualization_6_', dsN, '_', nCl, '.mat'];
fileForVisualization7Layer = [root, 'statistics/fileForVisualization_5_', dsN, '_', nCl, '.mat'];
fileForVisualization8Layer = [root, 'statistics/fileForVisualization_6_', dsN, '_', nCl, '.mat'];

% folders with visualized vocabulary
diskFolder = [commonRoot, 'Visualized vocabulary/'];
folderFor3Layer = [diskFolder, dataSetNames{dataSetNumber},'/3_layer/'];
folderFor4Layer = [diskFolder, dataSetNames{dataSetNumber},'/4_layer/'];
folderFor5Layer = [diskFolder, dataSetNames{dataSetNumber},'/5_layer/'];
folderFor6Layer = [diskFolder, dataSetNames{dataSetNumber},'/6_layer/'];
folderFor7Layer = [diskFolder, dataSetNames{dataSetNumber},'/7_layer/'];
folderFor8Layer = [diskFolder, dataSetNames{dataSetNumber},'/8_layer/'];

elPath2 = [depthPath, '_layer2'];
elPath3 = [depthPath, '_layer3'];
elPath4 = [depthPath, '_layer4'];
elPath5 = [depthPath, '_layer5'];
elPath6 = [depthPath, '_layer6'];
elPath7 = [depthPath, '_layer7'];
elPath8 = [depthPath, '_layer8'];

displ3 = 6;
displ5 = 18;
displ7 = 52;
isFIG = false;

%% define all parameters here

fileListPrecomputed = false;
is_subset = false; % whether we shall use all files for learning

% define the subset length
if dataSetNumber == 1 || dataSetNumber == 3
    subset_len = 10; % how much shall we use for training
    subsetPercent = 1.0; % not used
elseif dataSetNumber == 2
    subset_len = 1;
    subsetPercent = 0.8; % what percent from each folder to use
end

[dxKernel, combs, largestLine, sigma, sigmaKernelSize, isErrosion, discRadius, is_guided, r_guided, eps, ...
    is_mask_extended, maxExtThresh1, maxExtThresh2] = loadFilteringParameters(dataSetNumber);


%% Define what we have to learn

is_first_layer = false; % computes cluster centres, thresh and clusterSizes
is_third_layer = false;  % learns the first layer
is_4th_layer = false;
is_5th_layer = false; 
is_6th_layer = false;
is_7th_layer = true; 
is_8th_layer = false; 

is_statistics3_collected = false;
is_statistics4_collected = false;
is_statistics5_collected = false;
is_statistics6_collected = false;
is_statistics7_collected = false;
is_statistics8_collected = false;

is_statistics3_sieved = false;
is_statistics4_sieved = false;
is_statistics5_sieved = false;
is_statistics6_sieved = false;
is_statistics7_sieved = false;
is_statistics8_sieved = false;

is_statistics3_aggregated = false;
is_statistics4_aggregated = false;
is_statistics5_aggregated = false;
is_statistics6_aggregated = false;
is_statistics7_aggregated = false;
is_statistics8_aggregated = false;

is_partSelectionDone3 = false;
is_partSelectionDone4 = false;
is_partSelectionDone5 = false;
is_partSelectionDone6 = false;
is_partSelectionDone7 = false;
is_partSelectionDone8 = false;

% do we want to visualize the layer
visualizeLayer_3 = true;
visualizeLayer_4 = true;
visualizeLayer_5 = true;
visualizeLayer_6 = true;
visualizeLayer_7 = true;
visualizeLayer_8 = true;

%% here we define the first layer quantiles (bins of the first layer)

[quantilesFirst] = defineFirstLayerQuantiles(nClusters, dataSetNumber, is_guided);

%% creating a filelist here

if dataSetNumber == 1 || dataSetNumber == 3
    [list_depth, lenF] = extractFileList(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subset_len);
    list_mask = [];
elseif dataSetNumber == 2
    [list_depth, list_mask, ~, lenF] = extractFileListWashington(fileListPrecomputed, depthPath, depthPathDefault, is_subset, subsetPercent);
end

% % %  downsampling (if necessary)
% if dataSetNumber == 1 || dataSetNumber == 3
%     is_downsampling = true;
%     dowsample_rate = 1.3;
% elseif dataSetNumber == 2    
%     is_downsampling = true;
%     dowsample_rate = 2.5;
% end
% 
% 
% upsampleImages(list_depth, list_mask, lenF, is_downsampling, dowsample_rate, dataSetNumber); % to be done only once!


is_downsampling = false;
dowsample_rate = 1.0;

% fieldSize = [17, 5, 71];  % x, y and z directions
% [coverageOverall, areaOverall] = ProjectAllStatistics3(statistics3LayerSieved, statistics3LayerAggregated, fieldSize, ...
%                                                        list_depth, list_mask, lenF, dataSetNumber);
% ratio = coverageOverall/areaOverall;
% save(coverageOverall3Layer,  'coverageOverall');



%% learn the first layer

if is_first_layer  % here we learn the first layer parameters
    disp('Calibration of the first layer parameters...');

    [cluster1Centres, cluster1Bounds, thresh] = learnFirstLayer(list_depth, list_mask, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius, ...
                      nClusters, quantilesFirst, dataSetNumber, is_guided, r_guided, eps, ...
                      is_mask_extended, maxExtThresh1, maxExtThresh2);

    save(statistics1Layer,  'cluster1Centres', 'cluster1Bounds', 'thresh', 'nClusters', 'dataSetNumber');
end

% outFolder = '/home/vvk201/1LayerWash/';
% 
% applyFirstLayer(list_depth, list_mask, lenF, statistics1Layer, outFolder, sigma, sigmaKernelSize, ...
%                                 dxKernel, isErrosion, discRadius, is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2);



%% learn the third layer

if is_third_layer  % here we learn the third layer of the hierarchy
    
    disp('Learning of the third layer ...');
    % we have to perform the following procedures
    
    if ~is_first_layer
        % read the first layer
        load(statistics1Layer);
    end
    
    %load('displacements_Layers3_4.mat');
    load('settings/displacements_Layer3.mat'); % 'displacements'
    numDisps = size(displacements, 1);     % number of displacements
    fieldSize = [17, 7, 71];  % x, y and z directions
    quant = 0.06;
    borderOffset34 = 9;
    iterations = 400;
    depthStep = thresh/4;
    maxRelDepth = 127;
    
    
    elPath = elPath2;
    strFolderLen = length(depthPath);
    
    list_els = list_depth;
    for i = 1:lenF
        str = list_els{i};
        str = [elPath ,str(strFolderLen + 1:end)];
        list_els{i} = str;
    end
    

    if ~is_statistics3_collected
        
        [outputStatistics, outputCoords, curTS] = CollectStats_NextLayers(list_els, list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius, nClusters, ...
                                                        displacements, numDisps, borderOffset34, 3, depthStep, dataSetNumber, ...
                                                        is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, maxRelDepth);                                              
        save(statistics3Layer, 'outputStatistics', 'outputCoords');
    end
    
    if ~is_statistics3_sieved 
    
        load(statistics3Layer);
        numDisps = 2;
        thresh3Pair = 0.01 * lenF;
        [statistics, outputCoords, cluster3Depths, curTS] = sieveStatistics(outputStatistics, outputCoords, numDisps, thresh3Pair) 
        save(statistics3LayerSieved, 'statistics', 'cluster3Depths', 'outputCoords');      
    else
        
        load(statistics3LayerSieved);  %  variable statistics should be read
        curTS = size(statistics, 1);
    end
    
    if ~is_statistics3_aggregated % aggregate statistics
        
        sieve_thresh = 6;
        [X, frequencies, curTS, triples] = Aggregate_3Layer(statistics, n2Clusters, curTS, sieve_thresh);
        [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieve_thresh);
        outputCoords = outputCoords(ind, :); 
        
        save(statistics3LayerSieved, 'statistics', 'cluster3Depths', 'outputCoords');
        save(statistics3LayerAggregated, 'X' ,'frequencies', 'curTS', 'triples');
    end
    
    if ~is_partSelectionDone3
        if partSelectionMethod == 1
            
            lenSelected = 150;
            [triples3Out, coverageOut, n3Clusters] = PartSelectionFull(nClusters, n2Clusters, statistics3LayerSieved, statistics3LayerAggregated, ...
                    dataSetNumber, fieldSize, list_depth, lenF, meargeThresh3, iterations, 3, [], lenSelected);
            save(parts3Layer, 'triples3Out', 'coverageOut', 'n3Clusters');

        elseif partSelectionMethod == 2

            load(statistics3LayerSieved);    
            load(statistics3LayerAggregated);

            [triples3Out] = Optimization_layer3(X, frequencies, triples, nClusters, alphaParam(3), betaParam(3), gammaParam(3), subset_len); % this is the main procedure
            save('statistics/layer3.mat', 'triples3Out');  
        end
    end
    
    % store the visualization of the vocabulary to the folder
    
    if visualizeLayer_3
        
        load(parts3Layer);                 % 'triples3Out', 'coverageOut',    'n3Clusters'
        load(statistics3LayerSieved);      %  'statistics', 'cluster3Depths', 'outputCoords'
        
        fieldSize = [17, 17, 71];
        % load(fileForVisualization3Layer);

        % prepare data for visualization (convert to the right format)
        triple3OutDepth = store3Layer(triples3Out, cluster3Depths, n3Clusters, nClusters, partSelectionMethod);
        save(fileForVisualization3Layer, 'triple3OutDepth');

        [is_ok] = layer_N_demonstrator(3, [], [], [], [], [], triple3OutDepth, displ3, 0, 0, ...
                                            nClusters, n3Clusters, str_folder, fieldSize, depthStep, cluster1Centres, isFIG);
    end
    
end




%% learn the 4th layer

if is_4th_layer

    disp('Learning of the 4th layer ...');

    if ~is_first_layer % read the first layer
        load(statistics1Layer);
    end

    if ~is_third_layer % read the third layer
        load(parts3Layer); % 'triples3Out', 'coverageOut', 'n3Clusters');  
    end

    load('displacements_Layer4.mat');   % 'displacements'
    numDisps = size(displacements, 1); 
    borderOffset34 = 9;
    fieldSize = [17, 17, 71];
    depthStep = thresh/4;
    quant = 0.06;
    maxDist = 2;
    iterations = 300;
    maxRelDepth = 127;
    
    elPath = elPath3;
    strFolderLen = length(depthPath);
    
    list_els = list_depth;
    for i = 1:lenF
        str = list_els{i};
        str = [elPath ,str(strFolderLen + 1:end)];
        list_els{i} = str;
    end
    
    if ~is_statistics4_collected
        [outputStatistics, outputCoords, curTS] = CollectStats_NextLayers(list_els, list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius, nClusters, ...
                                displacements, numDisps, borderOffset34, 4, depthStep, dataSetNumber, ...
                                is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, maxRelDepth);
        save(statistics4Layer, 'outputStatistics', 'outputCoords');
    end

    if ~is_statistics4_sieved 
    
        load(statistics4Layer);
        numDisps = 2;
        thresh4Pair = 0.01 * lenF;
        [statistics, outputCoords, cluster4Depths, curTS] = sieveStatistics(outputStatistics, outputCoords, numDisps, thresh4Pair) 
        save(statistics4LayerSieved, 'statistics', 'cluster4Depths', 'outputCoords'); 
        
    else
        load(statistics4LayerSieved);  %  variable statistics should be read
        curTS = size(statistics, 1);
    end
    
    if ~is_statistics4_aggregated % aggregate statistics
        
        sieve_thresh = 3;
        [X, frequencies, curTS, triples] = Aggregate_3Layer(statistics, n3Clusters, curTS, sieve_thresh);
        
        % sieve 'statistics' and 'outputCoords' once again
        [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieve_thresh);
        outputCoords = outputCoords(ind, :);  
        
        save(statistics4LayerSieved, 'statistics', 'cluster4Depths', 'outputCoords');
        save(statistics4LayerAggregated, 'X' ,'frequencies', 'curTS', 'triples');
    end
    
    if ~is_partSelectionDone4
        if partSelectionMethod == 1 % libhop-based style 
            
            lenSelected = 250;
            [triples4Out, coverageOut, n4Clusters] = PartSelectionFull(nClusters, n2Clusters, statistics4LayerSieved, statistics4LayerAggregated, ...
                                     dataSetNumber, fieldSize, list_depth, lenF, meargeThresh4, ...
                                     iterations, 4, fileForVisualization3Layer, lenSelected);
                                 
        elseif partSelectionMethod == 2
            
            [triples4Out] = Optimization_layer4(X, frequencies, triples3Out, triples4, nClusters, n3Clusters, maxDist, subset_len, alphaParam(4), betaParam(4), gammaParam(4));
            n4Clusters = size(triples4Out, 1);
            save('statistics/layer4Wash1.mat', 'triples4Out');
        end
    end
    
    if visualizeLayer_4
        
        load(parts4Layer);                 % 'triples4Out', 'coverageOut',    'n4Clusters'
        load(statistics4LayerSieved);      %  'statistics', 'cluster4Depths', 'outputCoords'
        
        % prepare data for visualization (convert to the right format)
        [triple4OutDepth] = store4Layer(triples4Out, cluster4Depths, n4Clusters, nClusters, partSelectionMethod);
        
        load(fileForVisualization3Layer);  % 'triple3OutDepth'
        save(fileForVisualization4Layer, 'triple3OutDepth', 'triple4OutDepth');
        
        fieldSize = [17, 17, 71];
        [is_ok] = layer_N_demonstrator(4, [], [], [], [], triple4OutDepth, triple3OutDepth, displ3, 0, 0, ...
                                            nClusters, n4Clusters, str_folder, fieldSize, depthStep, cluster1Centres, isFIG);
    end

end   




%% learn the 5th layer

if is_5th_layer

    disp('Learning of the 5th layer ...');
    
    % read all the previous layers

    load(statistics1Layer);  % 'thresh', etc
    load(parts3Layer); % 'triples3Out', 'coverageOut', 'n3Clusters');  
    load(parts4Layer); % 'triples4Out', 'coverageOut', 'n4Clusters');
    

    load('displacements_Layer5.mat');   % 'displacements'
    numDisps = size(displacements, 1); 
    displacement56 = 26;
    fieldSize = [53, 17, 251];
    depthStep = thresh/4;
    quant = 0.04;
    maxDist = 2;
    iterations = 400;
    maxRelDepth = 500;
    
    elPath = elPath4;
    strFolderLen = length(depthPath);
    
    list_els = list_depth;
    for i = 1:lenF
        str = list_els{i};
        str = [elPath ,str(strFolderLen + 1:end)];
        list_els{i} = str;
    end
    
    if ~is_statistics5_collected
        [outputStatistics, outputCoords, curTS] = CollectStats_NextLayers(list_els, list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius, nClusters, ...
                                displacements, numDisps, displacement56, 5, depthStep, dataSetNumber, ...
                                is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, maxRelDepth);

        save(statistics5Layer, 'outputStatistics', 'outputCoords');
    end

    if ~is_statistics5_sieved 
    
        load(statistics5Layer);
        numDisps = 2;
        thresh5Pair = 15;
        [statistics, outputCoords, cluster5Depths, curTS] = sieveStatistics(outputStatistics, outputCoords, numDisps, thresh5Pair) 
        save(statistics5LayerSieved, 'statistics', 'cluster5Depths', 'outputCoords');       
    else
        load(statistics5LayerSieved);  %  variable statistics should be read
        curTS = size(statistics, 1);
    end
    
    if ~is_statistics5_aggregated % aggregate statistics
        
        sieve_thresh = 3;
        [X, frequencies, curTS, triples] = Aggregate_3Layer(statistics, n4Clusters, curTS, sieve_thresh);
        
        % sieve 'statistics' and 'outputCoords' once again
        [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieve_thresh);
        outputCoords = outputCoords(ind, :);  
        
        save(statistics5LayerSieved, 'statistics', 'cluster5Depths', 'outputCoords');
        save(statistics5LayerAggregated, 'X' ,'frequencies', 'curTS', 'triples');
    end
    
    if ~is_partSelectionDone5
        if partSelectionMethod == 1 % libhop-based style 
            
            lenSelected = 240;
            [triples5Out, coverageOut, n5Clusters] = PartSelectionFull(nClusters, n2Clusters, statistics5LayerSieved, statistics5LayerAggregated, ...
                                     dataSetNumber, fieldSize, list_depth, lenF, meargeThresh5, ...
                                     iterations, 5, fileForVisualization4Layer, lenSelected);                    
            save(parts5Layer, 'triples5Out', 'coverageOut', 'n5Clusters');
            
        elseif partSelectionMethod == 2  % to change
            
            disp('Not implemented');
        end
    end
    
    if visualizeLayer_5
        
        load(parts5Layer);                 % 'triples5Out', 'coverageOut',    'n5Clusters'
        load(statistics5LayerSieved);      % 'statistics', 'cluster5Depths', 'outputCoords'
        

        % prepare data for visualization (convert to the right format)
        [triple5OutDepth] = store4Layer(triples5Out, cluster5Depths, n5Clusters, nClusters, partSelectionMethod);
        load(fileForVisualization4Layer);  % 'triple3OutDepth', 'triple4OutDepth'
        save(fileForVisualization5Layer, 'triple3OutDepth', 'triple4OutDepth', 'triple5OutDepth');
        save(parts5Layer, 'triples5Out', 'coverageOut', 'n5Clusters');
        
        fieldSize = [53, 53, 251];
                      
        [is_ok] = layer_N_demonstrator(5, [], [], [], triple5OutDepth, triple4OutDepth, triple3OutDepth, displ3, displ5, 0, ...
                                            nClusters, n5Clusters, str_folder, fieldSize, depthStep, cluster1Centres, isFIG);
    end
end





%% learn the 6th layer

if is_6th_layer

    disp('Learning of the 6th layer ...');
    
    % read all the previous layers

    load(statistics1Layer);  % 'thresh', etc  
    
    load(parts5Layer); % 'triples5Out', 'coverageOut', 'n5Clusters');

    load('displacements_Layer6.mat');   % 'displacements'
    numDisps = size(displacements, 1); 
    displacement56 = 26;
    fieldSize = [53, 53, 351];
    depthStep = thresh/4;
    quant = 0.04;
    maxDist = 2;
    iterations = 440;
    maxRelDepth = 500;
    
    elPath = elPath5;
    strFolderLen = length(depthPath);
    
    list_els = list_depth;
    for i = 1:lenF
        str = list_els{i};
        str = [elPath ,str(strFolderLen + 1:end)];
        list_els{i} = str;
    end
    
    if ~is_statistics6_collected
        [outputStatistics, outputCoords, curTS] = CollectStats_NextLayers(list_els, list_depth, list_mask, lenF, sigma, sigmaKernelSize,...
                                dxKernel, isErrosion, discRadius, nClusters, ...
                                displacements, numDisps, displacement56, 6, depthStep, dataSetNumber, ...
                                is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, maxRelDepth);

        save(statistics6Layer, 'outputStatistics', 'outputCoords');
    end

    if ~is_statistics6_sieved 
    
        load(statistics6Layer);
        numDisps = 2;
        thresh6Pair = 8;
        [statistics, outputCoords, cluster6Depths, curTS] = sieveStatistics(outputStatistics, outputCoords, numDisps, thresh6Pair) 
        save(statistics6LayerSieved, 'statistics', 'cluster6Depths', 'outputCoords');       
    else
        load(statistics6LayerSieved);  %  variable statistics should be read
        curTS = size(statistics, 1);
    end
    
    if ~is_statistics6_aggregated % aggregate statistics
        
        sieve_thresh = 2;
        [X, frequencies, curTS, triples] = Aggregate_3Layer(statistics, n5Clusters, curTS, sieve_thresh);
        
        % sieve 'statistics' and 'outputCoords' once again
        [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieve_thresh);
        outputCoords = outputCoords(ind, :);
        
        save(statistics6LayerSieved, 'statistics', 'cluster6Depths', 'outputCoords');
        save(statistics6LayerAggregated, 'X' ,'frequencies', 'curTS', 'triples');
    end
    
    if ~is_partSelectionDone6
        if partSelectionMethod == 1 % libhop-based style 
            
            lenSelected = 400;
            [triples6Out, coverageOut, n6Clusters] = PartSelectionFull(nClusters, n2Clusters, statistics6LayerSieved, statistics6LayerAggregated, ...
                                     dataSetNumber, fieldSize, list_depth, lenF, meargeThresh6, ...
                                     iterations, 6, fileForVisualization5Layer, lenSelected);                 
            save(parts6Layer, 'triples6Out', 'coverageOut', 'n6Clusters');
            
        elseif partSelectionMethod == 2  % to change
            
            disp('Not Implemented!');
        end
    end
    
    if visualizeLayer_6
        
        load(parts6Layer);                 % 'triples6Out', 'coverageOut',    'n6Clusters'
        load(statistics6LayerSieved);      % 'statistics', 'cluster6Depths', 'outputCoords'
        
        
        % prepare data for visualization (convert to the right format)
        [triple6OutDepth] = store4Layer(triples6Out, cluster6Depths, n6Clusters, nClusters, partSelectionMethod);
        
        load(fileForVisualization5Layer);  % 'triple3OutDepth', 'triple4OutDepth' 'triple5OutDepth'
        save(fileForVisualization6Layer, 'triple3OutDepth', 'triple4OutDepth', 'triple5OutDepth', 'triple6OutDepth');
        
        fieldSize = [53, 53, 301];
        
        [is_ok] = layer_N_demonstrator(6, [], [], triple6OutDepth, triple5OutDepth, triple4OutDepth, triple3OutDepth, displ3, displ5, 0, ...
                                            nClusters, n6Clusters, str_folder, fieldSize, depthStep, cluster1Centres, isFIG); 
    end

end 




%% learn the 7th layer

if is_7th_layer

    disp('Learning of the 6th layer ...');
    
    % read all the previous layers

    load(statistics1Layer);  % 'thresh', etc  
    
    load(parts6Layer); % 'triples5Out', 'coverageOut', 'n5Clusters');

    load('displacements_Layer7.mat');   % 'displacements'
    numDisps = size(displacements, 1); 
    borderOffset78 = 54;
    fieldSize = [160, 160, 251];
    depthStep = thresh/4;
    quant = 0.04;
    maxDist = 2;
    iterations = 500;
    maxRelDepth = 1000;
    
    elPath = elPath6;
    strFolderLen = length(depthPath);
    
    list_els = list_depth;
    for i = 1:lenF
        str = list_els{i};
        str = [elPath ,str(strFolderLen + 1:end)];
        list_els{i} = str;
    end
    
    if ~is_statistics7_collected
        [outputStatistics, outputCoords, curTS] = CollectStats_NextLayers(list_els, list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius, nClusters, ...
                                displacements, numDisps, borderOffset78, 7, depthStep, dataSetNumber, ...
                                is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, maxRelDepth);

        save(statistics7Layer, 'outputStatistics', 'outputCoords');
    end

    if ~is_statistics7_sieved 
    
        load(statistics7Layer);
        numDisps = 2;
        thresh7Pair = 0.0005 * lenF;
        [statistics, outputCoords, cluster7Depths, curTS] = sieveStatistics(outputStatistics, outputCoords, numDisps, thresh7Pair);
        save(statistics7LayerSieved, 'statistics', 'cluster7Depths', 'outputCoords');
        
    else
        load(statistics7LayerSieved);  %  variable statistics should be read
        curTS = size(statistics, 1);
    end
    
    if ~is_statistics7_aggregated % aggregate statistics
        
        sieve_thresh = 2;
        [X, frequencies, curTS, triples] = Aggregate_3Layer(statistics, n6Clusters, curTS, sieve_thresh);
        
        % sieve 'statistics' and 'outputCoords' once again
        [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieve_thresh);
        outputCoords = outputCoords(ind, :);
        
        save(statistics7LayerSieved, 'statistics', 'cluster7Depths', 'outputCoords');
        save(statistics7LayerAggregated, 'X' ,'frequencies', 'curTS', 'triples');
    end
    
    if ~is_partSelectionDone7
        if partSelectionMethod == 1 % libhop-based style 
            
            lenSelected = 600;
            [triples7Out, coverageOut, n7Clusters] = PartSelectionFull(nClusters, n2Clusters, statistics7LayerSieved, statistics7LayerAggregated, ...
                                     dataSetNumber, fieldSize, list_depth, lenF, meargeThresh7, ...
                                     iterations, 7, fileForVisualization6Layer, lenSelected);
            save(parts7Layer, 'triples7Out', 'coverageOut', 'n7Clusters');
             
        elseif partSelectionMethod == 2  % to change 
            
            disp('Not implemented!');
        end
    end
    
    if visualizeLayer_7
        
        load(parts7Layer);                 % 'triples6Out', 'coverageOut',    'n6Clusters'
        load(statistics7LayerSieved);      % 'statistics', 'cluster6Depths', 'outputCoords'
        
        
        % prepare data for visualization (convert to the right format)
        [triple7OutDepth] = store4Layer(triples7Out, cluster7Depths, n7Clusters, nClusters, partSelectionMethod);
        
        load(fileForVisualization6Layer);  % 'triple3OutDepth', 'triple4OutDepth' 'triple5OutDepth'
        save(fileForVisualization7Layer, 'triple3OutDepth', 'triple4OutDepth', 'triple5OutDepth', 'triple6OutDepth', 'triple7OutDepth');
        
        fieldSize = [160, 160, 700];
                            
        [is_ok] = layer_N_demonstrator(7, [], triple7OutDepth, triple6OutDepth, triple5OutDepth, triple4OutDepth, triple3OutDepth, displ3, displ5, displ7, ...
                                    nClusters, n6Clusters, str_folder, fieldSize, depthStep, cluster1Centres, isFIG);                            
    end

end  



% matlabpool close 

