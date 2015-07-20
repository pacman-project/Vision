% This is a script to learn the 3D compositional hierarchical shape
% vocabulary

% dataSetNumber 
% for Aim@Shape dataSetNumber = 1;list_input
% for Washington dataSetNumber = 2;
% for Vladislav's_STD dataSetNumber = 3;
% This is the version LIBHOP3D_DG - with diffrential geometry


% inputDataType = 1; % depth images
% inputDataType = 2; % meshes
% inputDataType = 3; % point clouds

% partSelectionMethod
% for MDL-based for depth images:              partSelectionMethod = 1 
% for MDL-based for meshes and point clouds:   partSelectionMethod = 2


function [nNClusters] = learnHierarchy()

% define folders configuration
commonRoot = 'D:/';
root = [commonRoot, 'LibHoP3D_DG/'];
addPaths(root);

dataSetNumber = 5;
dataSetNames = getDataSetNames();

if dataSetNumber < 5
    inputDataType = 1;      % depth images
    partSelectionMethod = 1;
    
elseif dataSetNumber == 5
    inputDataType = 2;  % meshes
    partSelectionMethod = 2;
    
else
    inputDataType = 3;  % point clouds
    partSelectionMethod = 3;
end


nClusters = 7;
is_GPU_USED = false;
n2Clusters = nClusters^2;
if is_GPU_USED
    g = gpuDevice(1);
    reset(g);
end

is_multyScale = true;

if partSelectionMethod == 1 % for depth images only
    
     meargeThresh = defineMeargeThreshes(1, dataSetNumber);
     [iterations, lenSelected, fieldSize, maxRelDepth, quant, partsCoverArea, numSimilar, threshPair, sieve_thresh] = loadPartSelectionParameters(dataSetNumber);
    
elseif partSelectionMethod == 2  % for meshes
    
     meargeThresh = defineMeargeThreshes(1, dataSetNumber);
     [iterations, lenSelected, fieldSize, maxRelDepth, quant, partsCoverArea, numSimilar, threshPair, sieve_thresh] = loadPartSelectionParameters(dataSetNumber);
    
end

[inputPath] = getPathToData(dataSetNumber, commonRoot);
[scales, lineAdders]  = getScales(dataSetNumber, is_multyScale); 
computeAllscales(scales, inputPath{1,1}, lineAdders, dataSetNumber, inputDataType); % convert all input images to different scales



%% define output file names

dsN = num2str(dataSetNumber);
nCl = num2str(nClusters);

[partsLayer, fileForVisualizationLayer, partsLayerLoc, partsLayerEnt, partsSpecialSelected, partsLayerAll, calibrationFile] = getStandardFilePaths(root, dsN, nCl);

vocabulary1Layer = [root, 'statistics/statistics_1_', dsN, '_', nCl, '.mat'];
diskFolderVisVoc = [commonRoot, 'Visualized vocabulary/'];
diskFolderStatMap = [commonRoot, 'Visualized statistics/'];

for i = 3:8
    statisticalMapLayer{i}   = [root, 'statistics/statisticalMap_', num2str(i), '_', dsN, '_', nCl, '.mat'];
    folderForLayer{i}        = [diskFolderVisVoc,  dataSetNames{dataSetNumber},'/', num2str(i), '_layer/'];
    folderForLayerStatMap{i} = [diskFolderStatMap, dataSetNames{dataSetNumber},'/', num2str(i), '_layer/'];
end


%%
displ{3} = 6;
displ{5} = 18;
displ{7} = 52;
isFIG = false;

if ~exist(vocabulary1Layer, 'file');
    [cluster1Centres, cluster1Bounds, thresh] = defineFirstLayerQuantiles(nClusters, dataSetNumber);
    save(vocabulary1Layer, 'cluster1Centres', 'cluster1Bounds', 'thresh');
else
    load(vocabulary1Layer); % cluster1Centres, cluster1Bounds, thresh
end

% [ is_ok ] = VisualizeCornerParts(cluster1Centres); 
% a = 2;

%% define all parameters here
if inputDataType == 1   % depth images
%     [filtOptions, zeroThresh] = loadFilteringParameters(dataSetNumber);
%     options = GetOptions(filtOptions);
end

[learningElType, learningElRadius, ~, ~, offsetsConventional] = loadLearningInferenceStructElement(dataSetNumber);


%% Define what we have to learn
%                                    {1,2,3,4,5,6,7,8}
is_layer =                           {1,1,1,1,1,1,1,1};
is_statisticalMap =                  {0,0,0,0,0,0,0,0};
is_statistics_collection =           {0,0,1,1,1,1,0,0};

%                                    {1,2,3,4,5,6,7,8}
is_statistics_sieve_aggregate_Weak = {0,0,1,1,1,1,1,1};
is_localization =                    {0,0,0,0,0,0,0,0};
is_Entropy =                         {0,0,0,0,0,0,0,0};
is_partSelectionSpecialNeeded =      {0,0,0,0,0,0,0,0};

%                                    {1,2,3,4,5,6,7,8}
is_statistics_sieve_aggregate_Strong={0,0,1,1,1,0,0,0};
is_partSelectionNeeded =             {0,0,1,0,0,0,0,0};
combinePartSelection =               {0,0,1,1,1,1,0,0};
visualizeLayer =                     {0,0,1,1,1,1,1,1};
is_inferenceNeeded =                 {0,0,1,1,1,1,1,1};
isDownsampling  =                    {0,0,0,0,0,0,0,0};

%% downsampling of input images

is_downsampling = false;
if is_downsampling && inputDataType == 1
    dowsample_rate = 1.0;
    is_subsetDS = false; % whether we shall use all files for learning
    subsetPercentDS = 1.0;
    UpsampleAllImages(dataSetNumber, inputPath{1,1}, dowsample_rate, is_subsetDS, subsetPercentDS);
end


%% adding zero boundaries

% boundSize = 10;
% parfor i = 1:lenF
%     I = imread(list_input{i});
%     I = I(:,:,1);
%     I = addZeroBoundaries(I, boundSize);
%     imwrite(I, list_input{i}, 'png');
% end



%% inference of the first layer parts in the data

nNClusters{1} = nClusters;
nNClusters{2} = n2Clusters;

for i = 3:8
    tripleOutDepth{i}   = [];
    triplesCurOut{i}    = [];
    partsOutSpecial{i}   = [];
    nNClustersSpecial{i} = [];
    nNClusters{i} = 0;
end


depthStep = 0.0005;

% maxDist = 2;
% borderOffset34 = 9;
% thresh = 4;

%% calibrate the dataset using an image of the sphere-like objects

if inputDataType == 1  % depth images
    % perform calibration of the dataSet
    
    if exist(calibrationFile, 'file');
        zScale = load(calibrationFile);
        zScale = zScale.zScale;
    else
        % create this file
        zScale = calibrateDataSet(dataSetNumber, inputPath{1,1});
        save(calibrationFile, 'zScale');
    end
end

is_subset = false;
subsetPercent = 1.0;


%% learn the FIRST LEVEL PARTS

for layerID = 2:6
    
 if is_layer{layerID}
    
    nPrevClusters = nNClusters{layerID-1};
 
    if is_inferenceNeeded{layerID - 1} 
        infArray = zeros(1, 8);
        infArray(layerID-1) = 1; % inference of the previous layer
        LayersInference(infArray, dataSetNumber, inputDataType, nClusters, is_multyScale);
    end
    
    if layerID == 2
        continue;  % nothing to learn at this stage
    elseif layerID == 3
        % learn the third layer
        
    elseif layerID > 3
        load(partsLayerAll{layerID - 1});              % 'triplesCurOut', 'nNClusters', 'partsEntropy';
        load(fileForVisualizationLayer{layerID - 1});  % 'tripleOutDepth'
        clear('partsEntropy');
    end
              
%     if (layerID == 5 && isDownsampling{5}) || (layerID == 7 && isDownsampling{7}) % downsamplind required
%         
%         % take ALL depth images
%         is_subsetTemp = false;
%         subsetPercentTemp = 1;    
%         [list_inputAll, list_maskAll, ~, lenF] = extractFileListGeneral(inputPath{layerID-1, 1}, is_subsetTemp, subsetPercentTemp, dataSetNumber);
%         % make an elament list of ALL element images
%         [list_els] = makeElList(list_inputAll, inputPath{layerID-1, 1}, elPath{layerID-1, 2});
%         
%          % make depth images of the same resolution as downsampled marks images
%          downsampleImages(list_inputAll, inputPath{layerID-1, 1}, list_els, list_maskAll, dataSetNumber, inputPath{layerID-1, 2});  % downsamples depth images according to el images
%     end
%     
%     if ~strcmp(inputPath{layerID, 1}, inputPath{1,1})
%          % after this we have to take a new (downsampled) depth list  
%         [list_input, list_mask, ~, lenF] = extractFileListGeneral(inputPath{layerID, 1}, is_subset, subsetPercent, dataSetNumber);
%     end
    
%     % recompute list_els of elements from the previous layer
%     [list_els] = makeElList(list_input, inputPath{layerID, 1}, elPath{layerID-1, 2});
      
    str = ['Learning of the layer ', num2str(layerID), '...'];
    disp(str);    
 
    %% build statistical maps
    if is_statisticalMap{layerID}
        
        filterThresh = 20;
        [stats5D, sumSamples] = buildStatMap_NextLayers(list_els, list_input, list_mask, lenF, sigma, sigmaKernelSize, dxKernel, isErrosion, discRadius,...
                                                        nPrevClusters, depthStep, dataSetNumber, ...
                                                        is_guided, r_guided, eps,  is_mask_extended, maxExtThresh1, maxExtThresh2, fieldSize{layerID}, filterThresh);
                                                    
        save(statisticalMapLayer{layerID}, 'stats5D', 'sumSamples');
                                        
        % after this we may be interested to visualize the statistical maps
        visualizeStatisticalMaps( layerID - 1, tripleOutDepth, displ, ...
                                            nClusters, nPrevClusters, folderForLayerStatMap{layerID}, fieldSize{layerID}, depthStep, cluster1Centres, isFIG, stats5D);
    end
    
    
    %% collect statistics

    if is_statistics_collection{layerID} 

        for i = 1:length(scales)
            
            % get path to the input data at this scale
            str = inputPath{layerID, 1};
            input_path = getPathScale(str, lineAdders{i});
            elPath = getElPath(input_path, layerID-1);
            
            % get a list of elements and depth images
            [list_input, list_mask, ~, lenF] = extractFileListGeneral(input_path, is_subset, subsetPercent, dataSetNumber);
            [list_els] = makeElList(list_input, input_path, elPath);
            [statisticsLayer, ~, ~, ~, ~] = GetStatisticsFiles(root, dsN, nCl, lineAdders{i}, layerID);
            if inputDataType == 1
                distp('TODO');
            elseif inputDataType == 2
                [outputStatistics, outputCoords, outputFrames, curTS] = CollectStats_NextLayersMesh_PC(list_els, list_input, lenF, nPrevClusters, layerID, ...
                                                                        depthStep, dataSetNumber, maxRelDepth{3}, cluster1Bounds, nClusters, offsetsConventional);
            end
            
            save(statisticsLayer, 'outputStatistics', 'outputCoords', 'outputFrames', 'curTS', '-v7.3');
            save('Temp/depth_files.mat', 'list_input', 'list_mask', 'list_els');
        end
    end
    
    load('Temp/depth_files.mat');
    
    %% Weak sieve and aggregation of the statistics

    if is_statistics_sieve_aggregate_Weak{layerID}
        
        weak_multiplier = 0.1;
        
        for i = 1:length(scales)
            
            % load statistics
            [statisticsLayer, statisticsLayerSieved_Weak, statisticsLayerAggregated_Weak, ~, ~] = GetStatisticsFiles(root, dsN, nCl, lineAdders{i}, layerID);
            load(statisticsLayer);

            % sieve
            numDisps = 2;
            pairThresh = weak_multiplier *curTS * threshPair{layerID};
            [statistics, outputCoords, outputScales, outputFrames, clusterCurDepths, curTS] = sieveStatistics(outputStatistics, outputCoords, outputScales, outputFrames, numDisps, ...
                                                            pairThresh, nPrevClusters, quant{layerID}, maxRelDepth{layerID}, is_GPU_USED);
            clear('outputStatistics');

            % aggregate
            sieveThresh = weak_multiplier *curTS * sieve_thresh{layerID};
            is_GPU_USED = false;
            [X, frequencies, triples, is_sparse] = Aggregate_3Layer(statistics, nPrevClusters, curTS, sieveThresh, is_GPU_USED);
            [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieveThresh, is_sparse, is_GPU_USED);
            outputCoords = outputCoords(ind, :);
            outputScales = outputScales(ind, :);
            outputFrames = outputFrames(ind, :);
            curTS = size(outputCoords, 1);

            % save the processed statistics

            save(statisticsLayerSieved_Weak, 'statistics', 'clusterCurDepths', 'outputCoords', 'outputScales', 'outputFrames', '-v7.3');
            save(statisticsLayerAggregated_Weak, 'X' ,'frequencies', 'curTS', 'triples', '-v7.3');
        end
    end
    
    
  %% compute additional part selection criteria  
  
    if is_Entropy{layerID} % part selection based on discrimination criteria
        
        [partEntropy] = computePartsEntropy(list_input, statisticsLayerSieved_Weak{layerID}, statisticsLayerAggregated_Weak{layerID}, nPrevClusters+1, dataSetNumber);
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
        
%         if ~exist('partEntropy', 'var')
%             load('Temp/entropy.mat');
%         end
        
        winSize = 6;
        [localizationOut] = computePartsLocalization(list_input, statisticsLayerSieved_Weak{layerID}, statisticsLayerAggregated_Weak{layerID}, nPrevClusters+1, winSize);
%         save('Temp/localization.mat', 'localizationOut', 'partEntropy');           
        load(statisticsLayerAggregated_Weak{layerID});
        
%     Next lines come for visualization
        [localizationOut, inds] = sort(localizationOut, 'ascend');
        X = X(inds,:);
        a = [localizationOut,X];

        inds = localizationOut >= 0.25;
        triplesOutLoc = X(inds,:);
        numSelectedLoc = size(triplesOutLoc, 1);
        
        partLocEntropyOut = 5*ones(numSelectedLoc, 1);
        save(partsLayerLoc{layerID}, 'triplesOutLoc', 'numSelectedLoc', 'partLocEntropyOut');
        
    end
    
    
    %% localization-discriminative-based part selection procedure  
    
    if is_partSelectionSpecialNeeded{layerID}
        if ~exist('clusterCurDepths')
            bb = load(statisticsLayerSieved_Weak{layerID});
            clusterCurDepths = bb.clusterCurDepths;
        end
        [partsOutSpecial{layerID}, nNClustersSpecial{layerID}, partEntropySpecial] = partSelectionSpecial(false, true, partsLayerEnt{layerID},...
            partsLayerLoc{layerID},  nClusters, fileForVisualizationLayer, clusterCurDepths, displ, ...
            layerID, cluster1Centres, depthStep, partsCoverArea{layerID});          
  
        save(partsSpecialSelected{layerID}, 'partsOutSpecial', 'nNClustersSpecial', 'partEntropySpecial');
    end
    

    %% Strong Sieve and Aggregation of the statistics
    
    if is_statistics_sieve_aggregate_Strong{layerID} 
        
        strongMultiplier = 1.0;
        
        for i = 1:length(scales)
            
            [statisticsLayer, ~, ~, statisticsLayerSieved, statisticsLayerAggregated] = GetStatisticsFiles(root, dsN, nCl, lineAdders{i}, layerID);
            load(statisticsLayer);

            % sieve
            pairThresh = strongMultiplier * curTS * threshPair{layerID};
            numDisps = 2;

            [statistics, outputCoords, outputScales, outputFrames, clusterCurDepths, curTS] = sieveStatistics(outputStatistics, outputCoords, outputScales, outputFrames, numDisps, ...
                                                pairThresh, nPrevClusters, quant{layerID}, maxRelDepth{layerID}, is_GPU_USED);
            clear('outputStatistics');

            % aggregate        
            sieveThresh = strongMultiplier * curTS * sieve_thresh{layerID};
            [X, frequencies, triples, is_sparse] = Aggregate_3Layer(statistics, nPrevClusters, curTS, sieveThresh, is_GPU_USED);
            [ind, statistics] = Sieve3afterAggregation(statistics, triples, sieveThresh, is_sparse, is_GPU_USED);
            outputCoords = outputCoords(ind, :);
            outputScales = outputScales(ind, :);
            outputFrames = outputFrames(ind, :);
            curTS = size(outputCoords, 1);

            save(statisticsLayerSieved, 'statistics', 'clusterCurDepths', 'outputCoords', 'outputScales', 'outputFrames', '-v7.3');
            save(statisticsLayerAggregated, 'X' ,'frequencies', 'curTS', 'triples', '-v7.3');
        end
        
    end 
    
    
    
 %% Coverage-based part selection procedure   
         
    if is_partSelectionNeeded{layerID}
        if partSelectionMethod == 1
            
            if inputDataType == 1  % depth images
                [triplesCurOut{layerID}, coverageOut, nNClusters{layerID}] = PartSelectionFull(nClusters, nPrevClusters, dataSetNumber, partsCoverArea{layerID}, ... 
                    iterations{layerID}, layerID, fileForVisualizationLayer, lenSelected{layerID}, cluster1Centres, numSimilar{layerID}, offsetsConventional, depthStep, ...
                    is_multyScale, scales, lineAdders, inputPath{layerID, 1}, is_subset, subsetPercent, root, dsN, nCl);
            elseif inputDataType == 2  % meshes
                
            end
            
            partEntropyMain = 5*ones(nNClusters{layerID},1);    
            % previously selected 
            save(partsLayer{layerID}, 'triplesCurOut', 'coverageOut', 'nNClusters', 'partEntropyMain'); 
            
            
%             % define parts entropy: TODO
%             
%             load(partsLayer{layerID});
%             inFolder = elPath{layerID-1, 1};
%             outFolder = elPath{layerID, 1};
%             
%             InferenceNext_simple(statistics, outputCoords, outputScales, outputFrames, curTS, inFolder, outFolder, list_els, ...
%                                 triplesCurOut{layerID}, nNClusters{1}^2, nNClusters{layerID});
        end
    end
    
    % store the visualization of the vocabulary to the folder
    
    if combinePartSelection{layerID}
        
%         try
%             load(partsLayer{layerID});                 % 'triplesCurOut', 'coverageOut', 'nNClusters' 'partEntropyMain'
%             load(partsSpecialSelected{layerID});       % 'partsOutSpecial', 'nNClustersSpecial', 'partEntropySpecial'
%         
%         % this is to compute part's entropy for triplesCurOut
%         
%             partsEntropy{layerID} = [partEntropyMain; partEntropySpecial];
%             triplesCurOut{layerID} = [triplesCurOut{layerID}; partsOutSpecial{layerID}];
%             nNClusters{layerID} = nNClusters{layerID} + nNClustersSpecial{layerID};
% 
%             save(partsLayerAll{layerID}, 'triplesCurOut', 'nNClusters', 'partsEntropy');
%         catch exception
%     
%             % if there on of the files does not exist
%             load(partsLayer{layerID});
%             partsEntropy = partEntropyMain;
%             save(partsLayerAll{layerID}, 'triplesCurOut', 'nNClusters', 'partsEntropy');
%         end
           
        try
            load(partsSpecialSelected{layerID});       % 'partsOutSpecial', 'nNClustersSpecial', 'partEntropySpecial'
        
            % this is to compute part's entropy for triplesCurOut
        
            partsEntropy{layerID} = [partEntropySpecial];
            triplesCurOut{layerID} = [partsOutSpecial{layerID}];
            nNClusters{layerID} = nNClustersSpecial{layerID};

            save(partsLayerAll{layerID}, 'triplesCurOut', 'nNClusters', 'partsEntropy');
        catch exception
    
            % if there on of the files does not exist
            load(partsLayer{layerID});
            partsEntropy = partEntropyMain;
            save(partsLayerAll{layerID}, 'triplesCurOut', 'nNClusters', 'partsEntropy');
        end
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

        [is_ok] = layer_N_demonstrator(layerID, tripleOutDepth, offsetsConventional, ...
                                            nClusters, nNClusters{layerID}, folderForLayer{layerID}, fieldSize{LI}, depthStep, cluster1Centres, isFIG);
    end
    
  end
end

end
