% This is a script to learn the entire 3D compositional hierarchy

% dataSetNumber 
% for Aim@Shape dataSetNumber = 1;
% for Washington dataSetNumber = 2;
% for Vladislav's_STD dataSetNumber = 3;
% grasping dataset dataSetNumber = 4;


% partSelectionMethod
% for LibHob-motivated partSelectionMethod = 1 
% for optimization-based partSelectionMethod = 2



function [nNClusters] = learnHierarchy()

dataSetNumber = 3;
dataSetNames{1} = 'Aim@Shape';
dataSetNames{2} = 'Washington';
dataSetNames{3} = 'Vladislav_STD';

partSelectionMethod = 1;
nClusters = 7;

is_GPU_USED = false;
n2Clusters = nClusters^2;

% here we initialize the parallel computing (including GPU acceleration)
% g = gpuDevice(1);
% reset(g);
% matlabpool open 10    % for parallel computing


% define folders configuration
commonRoot = 'D:/';
root = [commonRoot, 'LibHoP3D/'];
addPaths(root);

if partSelectionMethod == 1
    
     meargeThresh = defineMeargeThreshes(1, dataSetNumber);
     [iterations, lenSelected, fieldSize, maxRelDepth, quant, partsCoverArea, numSimilar, threshPair, sieve_thresh] = loadPartSelectionParameters(dataSetNumber);
    
elseif partSelectionMethod == 2
    
    % parameters for the different layers
    alphaParam = [0, 0, 0, 0];
    betaParam =  [0, 0, 1.5, 1.0];
    gammaParam = [0, 0, 0.33, 1.0];
end


[depthPath, elPath] = getPathToData(dataSetNumber, commonRoot);


%% define output file names

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);

if partSelectionMethod == 1
    aL = '3';
end

[statisticsLayer, statisticsLayerSieved, statisticsLayerAggregated, statisticsLayerSieved_Weak, statisticsLayerAggregated_Weak, ...
    partsLayer, fileForVisualizationLayer, partsLayerLoc, partsLayerEnt, partsSpecialSelected, partsLayerAll] = getStandardFilePaths(root, dsN, nCl, aL);

% files with Statistical maps
statisticalMapLayer{3} = [root, 'statistics/statisticalMap_3_', dsN, '_', nCl, '.mat'];
statisticalMapLayer{4} = [root, 'statistics/statisticalMap_4_', dsN, '_', nCl, '.mat'];
statisticalMapLayer{5} = [root, 'statistics/statisticalMap_5_', dsN, '_', nCl, '.mat'];
statisticalMapLayer{6} = [root, 'statistics/statisticalMap_6_', dsN, '_', nCl, '.mat'];
statisticalMapLayer{7} = [root, 'statistics/statisticalMap_7_', dsN, '_', nCl, '.mat'];
statisticalMapLayer{8} = [root, 'statistics/statisticalMap_8_', dsN, '_', nCl, '.mat'];

% files for overall coverage
coverageOverallLayer{3} = [root, 'statistics/coverageOverall_3_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverallLayer{4} = [root, 'statistics/coverageOverall_4_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverallLayer{5} = [root, 'statistics/coverageOverall_5_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverallLayer{6} = [root, 'statistics/coverageOverall_6_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverallLayer{7} = [root, 'statistics/coverageOverall_7_', dsN, '_', nCl, '_a', aL, '.mat'];
coverageOverallLayer{8} = [root, 'statistics/coverageOverall_8_', dsN, '_', nCl, '_a', aL, '.mat'];

% folders with visualized vocabulary
diskFolder = [commonRoot, 'Visualized vocabulary/'];
folderForLayer{3} = [diskFolder, dataSetNames{dataSetNumber},'/3_layer/'];
folderForLayer{4} = [diskFolder, dataSetNames{dataSetNumber},'/4_layer/'];
folderForLayer{5} = [diskFolder, dataSetNames{dataSetNumber},'/5_layer/'];
folderForLayer{6} = [diskFolder, dataSetNames{dataSetNumber},'/6_layer/'];
folderForLayer{7} = [diskFolder, dataSetNames{dataSetNumber},'/7_layer/'];
folderForLayer{8} = [diskFolder, dataSetNames{dataSetNumber},'/8_layer/'];

% folders with visualized statistical maps
diskFolder = [commonRoot, 'Visualized statistics/'];
folderForLayerStatMap{3} = [diskFolder, dataSetNames{dataSetNumber},'/3_layer/'];
folderForLayerStatMap{4} = [diskFolder, dataSetNames{dataSetNumber},'/4_layer/'];
folderForLayerStatMap{5} = [diskFolder, dataSetNames{dataSetNumber},'/5_layer/'];
folderForLayerStatMap{6} = [diskFolder, dataSetNames{dataSetNumber},'/6_layer/'];
folderForLayerStatMap{7} = [diskFolder, dataSetNames{dataSetNumber},'/7_layer/'];
folderForLayerStatMap{8} = [diskFolder, dataSetNames{dataSetNumber},'/8_layer/'];


displ3 = 6;
displ5 = 18;
displ7 = 52;
isFIG = false;

%% define all parameters here


[dxKernel, dyKernelTop, dyKernelBottom, dxKernelBack, dxKernelForward, combs, largestLine, sigma, sigmaKernelSize, isErrosion, discRadius, is_guided, r_guided, eps, ...
    is_mask_extended, maxExtThresh1, maxExtThresh2] = loadFilteringParameters(dataSetNumber);
[learningElType, learningElRadius, ~, ~] = loadLearningInferenceStructElement(dataSetNumber);
[quantilesFirst] = defineFirstLayerQuantiles(nClusters, dataSetNumber, is_guided);


%% Define what we have to learn

is_layer{1} = 0; % computes cluster centres, thresh and clusterSizes
is_layer{3} = 0;  % learns the first layer
is_layer{4} = 1;
is_layer{5} = 0; 
is_layer{6} = 0;
is_layer{7} = 0; 
is_layer{8} = 0;

is_statisticalMap{3} = 0;
is_statisticalMap{4} = 0;
is_statisticalMap{5} = 0;
is_statisticalMap{6} = 0;
is_statisticalMap{7} = 0;
is_statisticalMap{8} = 0;

is_statistics_collection{3} = 1;
is_statistics_collection{4} = 1;
is_statistics_collection{5} = 0;
is_statistics_collection{6} = 0;
is_statistics_collection{7} = 0;
is_statistics_collection{8} = 0;

% Weak seive of statistics
is_statistics_sieve_aggregate_Weak{3} = 0;
is_statistics_sieve_aggregate_Weak{4} = 0;
is_statistics_sieve_aggregate_Weak{5} = 0;
is_statistics_sieve_aggregate_Weak{6} = 0;
is_statistics_sieve_aggregate_Weak{7} = 0;
is_statistics_sieve_aggregate_Weak{8} = 0;

% part selection according to special criteria
is_localization{3} = 0;
is_localization{4} = 0;
is_localization{5} = 0;
is_localization{6} = 0;
is_localization{7} = 0;
is_localization{8} = 0;

is_Entropy{3} = 0;
is_Entropy{4} = 0;
is_Entropy{5} = 0;
is_Entropy{6} = 0;
is_Entropy{7} = 0;
is_Entropy{8} = 0;

is_partSelectionSpecialNeeded{3} = 0;
is_partSelectionSpecialNeeded{4} = 0;
is_partSelectionSpecialNeeded{5} = 0;
is_partSelectionSpecialNeeded{6} = 0;
is_partSelectionSpecialNeeded{7} = 0;
is_partSelectionSpecialNeeded{8} = 0;

is_statistics_sieve_aggregate_Strong{3} = 1;
is_statistics_sieve_aggregate_Strong{4} = 1;
is_statistics_sieve_aggregate_Strong{5} = 0;
is_statistics_sieve_aggregate_Strong{6} = 0;
is_statistics_sieve_aggregate_Strong{7} = 0;
is_statistics_sieve_aggregate_Strong{8} = 0;

is_partSelectionNeeded{3} = 1;
is_partSelectionNeeded{4} = 1;
is_partSelectionNeeded{5} = 0;
is_partSelectionNeeded{6} = 0;
is_partSelectionNeeded{7} = 0;
is_partSelectionNeeded{8} = 0;

combinePartSelection{3} = 1;
combinePartSelection{4} = 1;
combinePartSelection{5} = 1;
combinePartSelection{6} = 1;
combinePartSelection{7} = 1;
combinePartSelection{8} = 1;

% do we want to visualize the layer
visualizeLayer{3} = 1;
visualizeLayer{4} = 1;
visualizeLayer{5} = 1;
visualizeLayer{6} = 1;
visualizeLayer{7} = 1;
visualizeLayer{8} = 1;

is_inferenceNeeded{2} = 1;
is_inferenceNeeded{3} = 1;
is_inferenceNeeded{4} = 1;
is_inferenceNeeded{5} = 1;
is_inferenceNeeded{6} = 0;
is_inferenceNeeded{7} = 0;
is_inferenceNeeded{8} = 0;

isDownsampling{5} = 0;
isDownsampling{7} = 0;

%% downsampling of ALL images

is_downsampling = false;
dowsample_rate = 1.0;

if is_downsampling  
    UpsampleAllImages(dataSetNumber, depthPath, dowsample_rate);
end


%% creating a filelist here

is_subset = false; % whether we shall use all files for learning

% define the subset length
if dataSetNumber == 1 || dataSetNumber == 3
    subset_len = 800; % how much shall we use for training
    subsetPercent = 1.0; % not used
elseif dataSetNumber == 2
    subset_len = 1;
    subsetPercent = 0.1; % what percent from each folder to use
end

if dataSetNumber == 1 || dataSetNumber == 3
    [list_depth, lenF] = extractFileList(depthPath{1}, is_subset, subset_len);
    list_mask = [];
elseif dataSetNumber == 2
    [list_depth, list_mask, ~, lenF] = extractFileListWashington(depthPath{1}, is_subset, subsetPercent);
end

%% adding zero boundaries


% boundSize = 10;
% parfor i = 1:lenF
%     I = imread(list_depth{i});
%     I = I(:,:,1);
%     I = addZeroBoundaries(I, boundSize);
%     imwrite(I, list_depth{i}, 'png');
% end



%% learn the first layer



if is_layer{1}  % here we learn the first layer parameters
    disp('Calibration of the first layer parameters...');

    [cluster1Centres, cluster1Bounds, thresh] = learnFirstLayer(list_depth, list_mask, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius, ...
                      nClusters, quantilesFirst, dataSetNumber, is_guided, r_guided, eps, ...
                      is_mask_extended, maxExtThresh1, maxExtThresh2);

    save(statisticsLayer{1},  'cluster1Centres', 'cluster1Bounds', 'thresh', 'nClusters', 'dataSetNumber');    
else
    load(statisticsLayer{1});
end


% outFolder = '/home/vvk201/1LayerWash/';
% 
% applyFirstLayer(list_depth, list_mask, lenF, statistics1Layer, outFolder, sigma, sigmaKernelSize, ...
%                                 dxKernel, isErrosion, discRadius, is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2);

nNClusters{1} = nClusters;
nNClusters{2} = n2Clusters;

for i = 3:8
    tripleOutDepth{i}   = [];
    triplesCurOut{i}    = [];
    partsOutSpecial{i}   = [];
    nNClustersSpecial{i} = [];
    nNClusters{i} = 0;
end

depthStep = thresh/4;
maxDist = 2;
borderOffset34 = 9;


%% learn the third layer

for layerID = 3:6
    
 if is_layer{layerID}  % here we learn the third layer of the hierarchy

    if mod(layerID, 2) == 1 
        load('settings/displacements_Layer3.mat');
    else
        load('settings/displacements_Layer4.mat');
    end
        
    if layerID ~= 3
        load(partsLayerAll{layerID - 1});              % 'triplesCurOut', 'nNClusters', 'partsEntropy';
        load(fileForVisualizationLayer{layerID - 1});  % 'tripleOutDepth'
        clear('partsEntropy');
    end
    
    nPrevClusters = nNClusters{layerID-1};
    
 
    if is_inferenceNeeded{layerID - 1}
        
        infArray = zeros(1, 8);
        infArray(layerID-1) = 1; % inference of the previous layer
        LayersInference(infArray, dataSetNumber, nClusters);
    end
    
    [list_els] = makeElList(list_depth, depthPath{layerID-1}, elPath{layerID-1});
              
    if layerID == 5 % downsamplind required

        if isDownsampling{5}  % make depth images of the same resolution as downsampled marks images
            downsampleImages(list_depth, depthPath, list_els, list_mask, dataSetNumber, [depthPath{4}, '_D1']);  % downsamples depth images according to el images
        end
    end
    
    if ~strcmp(depthPath{layerID}, depthPath{1})

        % create a file list again
        if dataSetNumber == 1 || dataSetNumber == 3
            [list_depth, lenF] = extractFileList(depthPath{layerID}, is_subset, subset_len);
            list_mask = [];
        elseif dataSetNumber == 2
            [list_depth, list_mask, ~, lenF] = extractFileListWashington(depthPath{layerID}, is_subset, subsetPercent);
        end
        
        % recompute list_els of elements from the previous layer
        [list_els] = makeElList(list_depth, depthPath{layerID-1}, elPath{layerID-1});
    end
    
    str = ['Learning of the layer ', num2str(layerID), '...'];
    disp(str);    
    numDisps = size(displacements, 1);     % number of displacements
    
 
    %% build statistical maps
    if is_statisticalMap{layerID}
        
        filterThresh = 20;
        [stats5D, sumSamples] = buildStatMap_NextLayers(list_els, list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius,...
                                                        nPrevClusters, depthStep, dataSetNumber, ...
                                                        is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, fieldSize{layerID}, filterThresh);
                                                    
        save(statisticalMapLayer{layerID}, 'stats5D', 'sumSamples');
                                        
        % after this we may be interested to visualize the statistical maps
        visualizeStatisticalMaps( layerID - 1, tripleOutDepth, displ3, displ5, displ7, ...
                                            nClusters, nPrevClusters, folderForLayerStatMap{layerID}, fieldSize{layerID}, depthStep, cluster1Centres, isFIG, stats5D);
    end
    
    
    %% collect statistics

    if is_statistics_collection{layerID}    
        
        [outputStatistics, outputCoords, curTS] = CollectStats_NextLayers(list_els, list_depth, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, ...
                                    discRadius, nPrevClusters, displacements, numDisps, layerID, depthStep, dataSetNumber, ...
                                    is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, maxRelDepth{layerID}, learningElType{layerID}, ...
                                    learningElRadius{layerID}, displ3, cluster1Bounds);
                                
        save(statisticsLayer{layerID}, 'outputStatistics', 'outputCoords', 'curTS', '-v7.3');
        save('Temp/depth_files.mat', 'list_depth', 'list_mask', 'list_els');
    end
    
    load('Temp/depth_files.mat');
    
    %% Weak sieve and aggregation of the statistics
    
    weak_multiplier = 0.1;
    
    if is_statistics_sieve_aggregate_Weak{layerID} 
        
        % sieve
        load(statisticsLayer{layerID});
        numDisps = 2;
        pairThresh = weak_multiplier *curTS * threshPair{layerID};
        [statistics, outputCoords, clusterCurDepths, curTS] = sieveStatistics(outputStatistics, outputCoords, numDisps, pairThresh, nPrevClusters, ...
                                                                quant{layerID}, maxRelDepth{layerID}, is_GPU_USED);

        % aggregate
        sieveThresh = weak_multiplier *curTS * sieve_thresh{layerID};
        [X, frequencies, triples, is_sparse] = Aggregate_3Layer(statistics, nPrevClusters, curTS, sieveThresh, is_GPU_USED);
        [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieveThresh, is_sparse, is_GPU_USED);
        outputCoords = outputCoords(ind, :);
        curTS = size(outputCoords, 1);
        
        save(statisticsLayerSieved_Weak{layerID}, 'statistics', 'clusterCurDepths', 'outputCoords', '-v7.3');
        save(statisticsLayerAggregated_Weak{layerID}, 'X' ,'frequencies', 'curTS', 'triples', '-v7.3');
        
        % if necessary we check how statistics is collected    
%         [coverageOverall, areaOverall] = ProjectAllStatistics(statisticsLayerSieved{layerID}, statisticsLayerAggregated{layerID}, partsCoverArea{layerID}, ...
%                                                                list_depth, list_mask, lenF, dataSetNumber);
%         ratio = coverageOverall/areaOverall;
    end
    
    
  %% compute additional part selection criteria  
  
    if is_Entropy{layerID} % part selection based on discrimination criteria
        
        [partEntropy] = computePartsEntropy(list_depth, statisticsLayerSieved_Weak{layerID}, statisticsLayerAggregated_Weak{layerID}, nPrevClusters+1, dataSetNumber);
        save('Temp/entropy.mat', 'partEntropy');
        
        load(statisticsLayerAggregated_Weak{layerID}); 
        
%     Next lines come for visualization
%         [partEntropy, inds] = sort(partEntropy, 'ascend');
%         X = X(inds,:);
%         a = [partEntropy, X];
        
        inds = partEntropy < 1.50;
        triplesOutEnt = X(inds,:);
        partEntropyOut = partEntropy(inds);
        numSelectedEnt = size(triplesOutEnt, 1);
        save(partsLayerEnt{layerID}, 'triplesOutEnt', 'numSelectedEnt', 'partEntropyOut');
    end
        
    if is_localization{layerID} % part selection based on localization criteria
        
        if ~exist('partEntropy', 'var')
            load('Temp/entropy.mat');
        end
        
        winSize = 6;
        [localizationOut] = computePartsLocalization(list_depth, statisticsLayerSieved_Weak{layerID}, statisticsLayerAggregated_Weak{layerID}, nPrevClusters+1, winSize);
        save('Temp/localization.mat', 'localizationOut', 'partEntropy');           
        load(statisticsLayerAggregated_Weak{layerID});
        
%     Next lines come for visualization
%         [localizationOut, inds] = sort(localizationOut, 'ascend');
%         X = X(inds,:);
%         a = [localizationOut,X];

        inds = localizationOut >= 0.40 & partEntropy < 1.6;
        triplesOutLoc = X(inds,:);
        numSelectedLoc = size(triplesOutLoc, 1);
        
        partLocEntropyOut = partEntropy(inds);
        save(partsLayerLoc{layerID}, 'triplesOutLoc', 'numSelectedLoc', 'partLocEntropyOut');
        
    end
    
    
    %% localization-discriminative-based part selection procedure  
    
    if is_partSelectionSpecialNeeded{layerID}
        
        [partsOutSpecial{layerID}, nNClustersSpecial{layerID}, partEntropySpecial] = partSelectionSpecial(true, true, partsLayerEnt{layerID}, partsLayerLoc{layerID},  ...
                nClusters, fileForVisualizationLayer, clusterCurDepths, displ3, displ5, displ7, ...
                layerID, cluster1Centres, depthStep, partsCoverArea{layerID});          
  
        save(partsSpecialSelected{layerID}, 'partsOutSpecial', 'nNClustersSpecial', 'partEntropySpecial');
    end
    

    %% Strong Sieve and Aggregation of the statistics
    
    if is_statistics_sieve_aggregate_Strong{layerID} 
        
        % sieve
        pairThresh = curTS * threshPair{layerID};
        load(statisticsLayer{layerID});
        numDisps = 2;
        [statistics, outputCoords, clusterCurDepths, curTS] = sieveStatistics(outputStatistics, outputCoords, numDisps, pairThresh, nPrevClusters, ...
                                                                            quant{layerID}, maxRelDepth{layerID}, is_GPU_USED);
        
        % aggregate        
        sieveThresh = curTS * sieve_thresh{layerID};
        [X, frequencies, triples, is_sparse] = Aggregate_3Layer(statistics, nPrevClusters, curTS, sieveThresh, is_GPU_USED);
        [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieveThresh, is_sparse, is_GPU_USED);
        outputCoords = outputCoords(ind, :); 
        
        save(statisticsLayerSieved{layerID}, 'statistics', 'clusterCurDepths', 'outputCoords', '-v7.3');
        save(statisticsLayerAggregated{layerID}, 'X' ,'frequencies', 'curTS', 'triples', '-v7.3');
        
        % if necessary we check how statistics is collected
        
%         [coverageOverall, areaOverall] = ProjectAllStatistics(statisticsLayerSieved{layerID}, statisticsLayerAggregated{layerID}, partsCoverArea{layerID}, ...
%                                                                list_depth, list_mask, lenF, dataSetNumber);
%         ratio = coverageOverall/areaOverall;
    end 
    
    
    
 %% Coverage-based part selection procedure   
         
    if is_partSelectionNeeded{layerID}
        if partSelectionMethod == 1
            
            [triplesCurOut{layerID}, coverageOut, nNClusters{layerID}] = PartSelectionFull(nClusters, nPrevClusters, statisticsLayerSieved{layerID},...
                    statisticsLayerAggregated{layerID}, dataSetNumber, partsCoverArea{layerID}, list_depth, lenF, iterations{layerID}, layerID, ...
                    fileForVisualizationLayer, lenSelected{layerID}, displ3, displ5, displ7, cluster1Centres, depthStep, numSimilar{layerID}, is_GPU_USED);
            
            % define parts entropy: TODO
            
            partEntropyMain = 5*ones(nNClusters{layerID},1);    
            % previously selected 
            save(partsLayer{layerID}, 'triplesCurOut', 'coverageOut', 'nNClusters', 'partEntropyMain'); 
        end
    end
    
    % store the visualization of the vocabulary to the folder
    
    if combinePartSelection{layerID}
        
        load(partsLayer{layerID});                 % 'triplesCurOut', 'coverageOut', 'nNClusters' 'partEntropyMain'
        partsEntropy = partEntropyMain; 
        
%         load(partsSpecialSelected{layerID});       % 'partsOutSpecial', 'nNClustersSpecial', 'partEntropySpecial'
%         
%         % this is to compute part's entropy for triplesCurOut
% %         TODO !!!!
%         partsEntropy{layerID} = [partEntropyMain; partEntropySpecial];
%         triplesCurOut{layerID} = [triplesCurOut{layerID}; partsOutSpecial{layerID}];
%         nNClusters{layerID} = nNClusters{layerID} + nNClustersSpecial{layerID};

        save(partsLayerAll{layerID}, 'triplesCurOut', 'nNClusters', 'partsEntropy');
    end
    
    
    %% visualization of parts of the layer
    
    if visualizeLayer{layerID}

        load(partsLayerAll{layerID});
        load(statisticsLayerSieved{layerID});      %  'statistics', 'clusterCurDepths', 'outputCoords'
        
        if mod(layerID,2) == 1
            LI = layerID + 1;  % to make a squared visualization field
        else
            LI = layerID;
        end
            
        % load(fileForVisualization3Layer);

        % prepare data for visualization (convert to the right format)
        if layerID == 3
            tripleOutDepth{layerID} = store3Layer(triplesCurOut{layerID}, clusterCurDepths, nNClusters{layerID}, nClusters, partSelectionMethod);
        else
            tripleOutDepth{layerID} = store4Layer(triplesCurOut{layerID}, clusterCurDepths, nNClusters{layerID}, nClusters, partSelectionMethod);
        end
        
        save(fileForVisualizationLayer{layerID}, 'tripleOutDepth');

        [is_ok] = layer_N_demonstrator(layerID, tripleOutDepth, displ3, displ5, displ7, ...
                                            nClusters, nNClusters{layerID}, folderForLayer{layerID}, fieldSize{LI}, depthStep, cluster1Centres, isFIG);
    end
    
  end
end

end


